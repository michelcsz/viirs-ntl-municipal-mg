# ================================================================
# PROVA DE CONCEITO, PIB MENSAL MODELADO PARA OS 23 MUNICIPIOS
# DO VALE DO MUCURI (MG) A PARTIR DE LUZES NOTURNAS VIIRS
# VERSAO COM FILTRO DE QUALIDADE, WINSORIZACAO E DIAGNOSTICOS
# ================================================================
# O que este script faz:
# 1) baixa as geometrias dos 23 municipios do Vale do Mucuri;
# 2) extrai estatisticas mensais de luminosidade VIIRS no Earth Engine;
# 3) baixa o PIB municipal anual oficial do IBGE/SIDRA (2012-2023);
# 4) estima um modelo anual simples: log(PIB) ~ log(luz anual) + FE municipal + tendencia;
# 5) construi uma proxy mensal benchmarked ao PIB anual oficial;
# 6) aplica penalizacao por baixa cobertura e winsorizacao da luz mensal;
# 7) gera diagnosticos automaticos, inclusive para Teofilo Otoni;
# 8) deixa um bloco opcional para Chow-Lin mensal.
#
# IMPORTANTE:
# - O resultado principal nao e PIB mensal oficial. E uma proxy/modelagem mensal
#   ancorada no PIB anual oficial do IBGE.
# - A serie ancorada vai de 2012-01 a 2023-12.
# - O bloco Chow-Lin e opcional e pode ser usado como upgrade metodologico.
#
# PRE-REQUISITOS:
# - Conta no Google Earth Engine (https://earthengine.google.com/)
# - Projeto GEE criado e autenticado via rgee
# - Python configurado com o ambiente rgee (ver documentacao do pacote rgee)
# ================================================================

# =========================================================
# 0) PARAMETROS GERAIS
# =========================================================

reticulate_python <- "caminho/para/seu/python/rgee"

# Substitua pelo ID do seu projeto no Google Earth Engine.
# O projeto e criado em https://console.cloud.google.com/
ee_project <- "seu-projeto-gee"

start_date <- "2012-01-01"
end_date_exclusive <- "2024-01-01"   # exclusivo, portanto vai ate 2023-12
sidra_years <- 2012:2023

save_outputs <- TRUE
run_chow_lin <- TRUE   # mude para TRUE quando quiser rodar a versao mensal econometrica

# Lista dos 23 municipios Vale do Mucuri
mucuri_names <- c(
  "Ataleia", "Catuji", "Franciscopolis", "Frei Gaspar", "Itaipe",
  "Ladainha", "Malacacheta", "Novo Oriente de Minas", "Ouro Verde de Minas",
  "Pavao", "Pote", "Setubinha", "Teofilo Otoni",
  "Aguas Formosas", "Bertopolis", "Carlos Chagas", "Crisolita",
  "Fronteira dos Vales", "Machacalis", "Nanuque", "Santa Helena de Minas",
  "Serra dos Aimores", "Umburatiba"
)

# =========================================================
# 1) PYTHON / EARTH ENGINE
# =========================================================
Sys.setenv(RETICULATE_PYTHON = reticulate_python)

library(reticulate)
reticulate::py_config()

cert_path <- reticulate::py_eval("__import__('certifi').where()")
Sys.setenv(
  SSL_CERT_FILE = cert_path,
  REQUESTS_CA_BUNDLE = cert_path,
  CURL_CA_BUNDLE = cert_path
)

# =========================================================
# 2) PACOTES
# =========================================================
pkgs <- c(
  "tidyverse", "lubridate", "sf", "geobr", "rgee", "sidrar",
  "fixest", "slider", "janitor", "stringi", "tempdisagg",
  "geojsonio", "zoo"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

invisible(lapply(pkgs, library, character.only = TRUE))

options(timeout = 600)
options(scipen = 999)

# =========================================================
# 3) FUNCOES AUXILIARES
# =========================================================
norm_txt <- function(x) {
  x %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    toupper() %>%
    stringr::str_squish()
}

winsorize_vec <- function(x, probs = c(0.01, 0.99)) {
  q <- stats::quantile(x, probs = probs, na.rm = TRUE, type = 7)
  x <- pmin(pmax(x, q[[1]]), q[[2]])
  x
}

ensure_rgee_sessioninfo <- function(email) {
  # Substitua pelo e-mail associado ao seu projeto no Google Earth Engine.
  ee_path <- path.expand("~/.config/earthengine")
  dir.create(ee_path, recursive = TRUE, showWarnings = FALSE)

  session_info <- data.frame(
    email = email,
    drive_cre = NA_character_,
    gcs_cre = NA_character_,
    stringsAsFactors = FALSE
  )

  write.table(
    session_info,
    file = file.path(ee_path, "rgee_sessioninfo.txt"),
    row.names = FALSE,
    quote = TRUE
  )
}

pick_col <- function(nms, patterns) {
  hit <- nms[stringr::str_detect(nms, paste(patterns, collapse = "|"))]
  if (length(hit) == 0) {
    stop("Nao encontrei a coluna esperada no retorno do SIDRA. Verifique names(pib_raw).")
  }
  hit[1]
}

estimate_monthly_td <- function(df_month, df_year) {
  df_month <- df_month %>% arrange(date)
  df_year  <- df_year  %>% arrange(year)

  fallback <- tibble(
    code_muni = first(df_month$code_muni),
    name_muni = first(df_month$name_muni),
    date = df_month$date,
    pib_month_td = df_month$pib_month_est,
    td_method = "fallback_proportional"
  )

  if (nrow(df_month) == 0 || nrow(df_year) == 0) {
    return(fallback)
  }

  x <- df_month$signal_smooth
  x <- ifelse(is.finite(x), x, NA_real_)
  x <- zoo::na.approx(x, na.rm = FALSE, rule = 2)

  if (all(is.na(x))) {
    return(fallback)
  }

  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  x <- pmax(x, 0) + 1e-6

  y_a <- ts(df_year$pib_annual, start = min(df_year$year), frequency = 1)
  x_m <- ts(
    x,
    start = c(min(df_month$year), min(df_month$month)),
    frequency = 12
  )

  fit_td <- tryCatch(
    tempdisagg::td(
      y_a ~ x_m,
      to = "monthly",
      method = "chow-lin-maxlog",
      conversion = "sum"
    ),
    error = function(e) NULL
  )

  if (is.null(fit_td)) {
    return(fallback)
  }

  y_m <- tryCatch(
    as.numeric(predict(fit_td)),
    error = function(e) NULL
  )

  if (is.null(y_m) || length(y_m) != nrow(df_month)) {
    return(fallback)
  }

  tibble(
    code_muni = first(df_month$code_muni),
    name_muni = first(df_month$name_muni),
    date = df_month$date,
    pib_month_td = y_m,
    td_method = "chow-lin-maxlog"
  )
}

# =========================================================
# 4) INICIALIZACAO DO EARTH ENGINE
# =========================================================
library(rgee)
ee$Initialize(project = ee_project)

# Substitua pelo e-mail associado ao seu projeto no Google Earth Engine.
ensure_rgee_sessioninfo(email = "seu-email@gmail.com")

# =========================================================
# 5) GEOMETRIAS MUNICIPAIS
# =========================================================
mg_muni <- geobr::read_municipality(
  code_muni = "MG",
  year = 2020,
  simplified = TRUE,
  showProgress = FALSE
)

mucuri_sf <- mg_muni %>%
  mutate(name_norm = norm_txt(name_muni)) %>%
  filter(name_norm %in% norm_txt(mucuri_names)) %>%
  st_transform(4326) %>%
  select(code_muni, name_muni, abbrev_state)

stopifnot(nrow(mucuri_sf) == 23)

muni_codes <- mucuri_sf$code_muni %>% as.integer()
mucuri_ee <- rgee::sf_as_ee(mucuri_sf, via = "getInfo")

# =========================================================
# 6) EARTH ENGINE / VIIRS MENSAL
# =========================================================
viirs_ic <- ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG")$
  filterDate(start_date, end_date_exclusive)

extract_viirs_month <- function(img) {
  date_ee <- ee$Date(img$get("system:time_start"))
  ym      <- date_ee$format("YYYY-MM")

  rad <- img$select("avg_rad")$max(0)$rename("rad")
  cf  <- img$select("cf_cvg")

  rad_masked <- rad$updateMask(cf$gt(0))
  stack_img  <- rad_masked$addBands(cf)

  reducer_all <- ee$Reducer$mean()$combine(
    reducer2 = ee$Reducer$sum(),
    sharedInputs = TRUE
  )

  out <- stack_img$reduceRegions(
    collection = mucuri_ee,
    reducer = reducer_all,
    scale = 500
  )$map(function(f) {
    f$set(
      "ym", ym,
      "date", date_ee$format("YYYY-MM-dd")
    )
  })

  out
}

monthly_fc <- ee$FeatureCollection(viirs_ic$map(extract_viirs_month)$flatten())
monthly_sf <- rgee::ee_as_sf(monthly_fc, via = "getInfo")

lights_monthly_raw <- monthly_sf %>%
  st_drop_geometry() %>%
  janitor::clean_names() %>%
  transmute(
    code_muni = as.integer(code_muni),
    name_muni = as.character(name_muni),
    date      = as.Date(date),
    ym        = as.character(ym),
    year      = lubridate::year(date),
    month     = lubridate::month(date),
    rad_mean  = as.numeric(rad_mean),
    rad_sum   = as.numeric(rad_sum),
    cf_cvg_mean = as.numeric(cf_cvg_mean),
    cf_cvg_sum  = as.numeric(cf_cvg_sum)
  ) %>%
  arrange(code_muni, date) %>%
  mutate(
    rad_mean = pmax(rad_mean, 0),
    rad_sum  = pmax(rad_sum, 0)
  )

# =========================================================
# 7) LIMPEZA / QUALIDADE DA LUZ MENSAL
# =========================================================
lights_monthly <- lights_monthly_raw %>%
  group_by(code_muni) %>%
  mutate(
    cov_ref = median(cf_cvg_mean[is.finite(cf_cvg_mean) & cf_cvg_mean > 0], na.rm = TRUE),
    cov_ref = ifelse(is.finite(cov_ref), cov_ref, NA_real_),
    quality_factor = case_when(
      is.na(cf_cvg_mean) ~ 0,
      is.na(cov_ref) ~ if_else(cf_cvg_mean > 0, 1, 0),
      cov_ref <= 0 ~ if_else(cf_cvg_mean > 0, 1, 0),
      TRUE ~ pmin(cf_cvg_mean / cov_ref, 1)
    ),
    rad_sum_win = winsorize_vec(rad_sum, probs = c(0.01, 0.99)),
    signal_raw = rad_sum_win * quality_factor
  ) %>%
  arrange(code_muni, date) %>%
  mutate(
    signal_smooth = slider::slide_dbl(
      signal_raw,
      ~ mean(.x, na.rm = TRUE),
      .before = 1,
      .after = 1,
      .complete = FALSE
    ),
    signal_smooth = if_else(is.finite(signal_smooth), pmax(signal_smooth, 0), NA_real_)
  ) %>%
  ungroup()

# =========================================================
# 8) PIB MUNICIPAL ANUAL DO SIDRA
# =========================================================
api_path <- paste0(
  "/t/5938/n6/", paste(muni_codes, collapse = ","),
  "/v/37/p/", paste(sidra_years, collapse = ","),
  "/d/v37%200"
)

pib_raw <- sidrar::get_sidra(api = api_path) %>%
  janitor::clean_names()

nms <- names(pib_raw)
col_code  <- pick_col(nms, c("^municipio_codigo$", "municipio_.*codigo", "munic.pio_.*codigo"))
col_name  <- pick_col(nms, c("^municipio$", "^munic.pio$"))
col_year  <- pick_col(nms, c("^ano$", "^periodo$"))
col_value <- pick_col(nms, c("^valor$"))

pib_annual <- pib_raw %>%
  transmute(
    code_muni = as.integer(.data[[col_code]]),
    name_muni = as.character(.data[[col_name]]),
    year      = as.integer(.data[[col_year]]),
    pib_annual = as.numeric(.data[[col_value]])
  ) %>%
  arrange(code_muni, year)

# =========================================================
# 9) AGREGACAO ANUAL DAS LUZES
# =========================================================
lights_annual <- lights_monthly %>%
  group_by(code_muni, name_muni, year) %>%
  summarise(
    light_year_sum_signal = sum(signal_smooth, na.rm = TRUE),
    light_year_sum_raw    = sum(rad_sum, na.rm = TRUE),
    light_year_mean_raw   = mean(rad_mean, na.rm = TRUE),
    cf_year_mean          = mean(cf_cvg_mean, na.rm = TRUE),
    .groups = "drop"
  )

panel_annual <- pib_annual %>%
  left_join(lights_annual, by = c("code_muni", "year"), suffix = c("_pib", "_light")) %>%
  mutate(
    name_muni = coalesce(name_muni_pib, name_muni_light),
    t = year - min(year) + 1,
    log_pib = log(pib_annual),
    log_light = log1p(light_year_sum_signal)
  ) %>%
  select(
    code_muni, name_muni, year, pib_annual,
    light_year_sum_signal, light_year_sum_raw, light_year_mean_raw,
    cf_year_mean, t, log_pib, log_light
  )

# =========================================================
# 10) MODELO ANUAL SIMPLES DE CALIBRACAO
# =========================================================
mod_annual <- fixest::feols(
  log_pib ~ log_light + t | code_muni,
  data = panel_annual
)

smearing <- mean(exp(residuals(mod_annual)), na.rm = TRUE)

panel_annual <- panel_annual %>%
  mutate(
    pib_hat_annual_in = exp(predict(mod_annual, newdata = panel_annual)) * smearing,
    abs_pct_error = abs(pib_hat_annual_in / pib_annual - 1) * 100
  )

annual_metrics <- panel_annual %>%
  summarise(
    rmse = sqrt(mean((pib_hat_annual_in - pib_annual)^2, na.rm = TRUE)),
    mape = mean(abs_pct_error, na.rm = TRUE),
    cor  = cor(pib_annual, pib_hat_annual_in, use = "complete.obs")
  )

print(summary(mod_annual))
print(annual_metrics)

# =========================================================
# 11) PROXY MENSAL BENCHMARKED AO PIB ANUAL
# =========================================================
monthly_shares <- lights_monthly %>%
  mutate(
    signal_smooth = if_else(is.finite(signal_smooth), signal_smooth, 0),
    signal_smooth = pmax(signal_smooth, 0)
  ) %>%
  group_by(code_muni, year) %>%
  mutate(
    yearly_signal_total = sum(signal_smooth, na.rm = TRUE),
    n_months_in_year = dplyr::n(),
    light_share_in_year = dplyr::if_else(
      yearly_signal_total > 0,
      signal_smooth / yearly_signal_total,
      1 / n_months_in_year
    )
  ) %>%
  ungroup()

pib_monthly <- monthly_shares %>%
  left_join(
    pib_annual %>% select(code_muni, year, pib_annual),
    by = c("code_muni", "year")
  ) %>%
  mutate(
    pib_month_est = pib_annual * light_share_in_year,
    source = "monthly_proxy_benchmarked_to_official_annual"
  ) %>%
  select(
    code_muni, name_muni, date, year, month,
    rad_mean, rad_sum, cf_cvg_mean, cf_cvg_sum,
    quality_factor, rad_sum_win, signal_raw, signal_smooth,
    light_share_in_year, pib_annual, pib_month_est, source
  ) %>%
  arrange(code_muni, date)

# =========================================================
# 12) CHECAGEM DE COERENCIA ANUAL
# =========================================================
check_benchmark <- pib_monthly %>%
  group_by(code_muni, name_muni, year) %>%
  summarise(
    pib_monthly_sum = sum(pib_month_est, na.rm = TRUE),
    pib_annual = first(pib_annual),
    diff = pib_monthly_sum - pib_annual,
    rel_diff = pib_monthly_sum / pib_annual - 1,
    .groups = "drop"
  )

print(summary(check_benchmark$diff))
print(summary(check_benchmark$rel_diff))

# =========================================================
# 13) DIAGNOSTICOS AUTOMATICOS, TEOFILO OTONI
# =========================================================
teofilo_diag <- pib_monthly %>%
  filter(norm_txt(name_muni) == "TEOFILO OTONI") %>%
  arrange(date) %>%
  mutate(
    z_rad_sum = (rad_sum - mean(rad_sum, na.rm = TRUE)) / sd(rad_sum, na.rm = TRUE),
    z_signal  = (signal_smooth - mean(signal_smooth, na.rm = TRUE)) / sd(signal_smooth, na.rm = TRUE)
  )

teofilo_outliers <- teofilo_diag %>%
  arrange(desc(abs(z_signal))) %>%
  select(
    date, year, month, rad_sum, cf_cvg_mean, quality_factor,
    signal_smooth, light_share_in_year, pib_annual, pib_month_est,
    z_rad_sum, z_signal
  ) %>%
  slice(1:12)

teofilo_seasonality <- teofilo_diag %>%
  group_by(month) %>%
  summarise(
    avg_light_share = mean(light_share_in_year, na.rm = TRUE),
    avg_signal = mean(signal_smooth, na.rm = TRUE),
    avg_pib_month = mean(pib_month_est, na.rm = TRUE),
    .groups = "drop"
  )

print(teofilo_outliers)
print(teofilo_seasonality)

# =========================================================
# 14) BLOCO OPCIONAL, CHOW-LIN MENSAL
# =========================================================
if (run_chow_lin) {
  chow_lin_monthly <- pib_monthly %>%
    group_split(code_muni) %>%
    purrr::map_dfr(function(d) {
      muni_id <- unique(d$code_muni)
      d_year  <- pib_annual %>% filter(code_muni == muni_id)
      estimate_monthly_td(d, d_year)
    })

  pib_monthly <- pib_monthly %>%
    left_join(
      chow_lin_monthly,
      by = c("code_muni", "name_muni", "date")
    )

  chow_lin_check <- pib_monthly %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(code_muni, name_muni, year) %>%
    summarise(
      pib_td_sum = sum(pib_month_td, na.rm = TRUE),
      pib_annual = first(pib_annual),
      diff_td = pib_td_sum - pib_annual,
      rel_diff_td = pib_td_sum / pib_annual - 1,
      .groups = "drop"
    )

  print(summary(chow_lin_check$diff_td))
  print(summary(chow_lin_check$rel_diff_td))
}

# =========================================================
# 15) SAIDAS
# =========================================================
if (save_outputs) {
  readr::write_csv(lights_monthly,      "lights_monthly_vale_mucuri_2012_2023_clean.csv")
  readr::write_csv(panel_annual,        "panel_annual_vale_mucuri_2012_2023.csv")
  readr::write_csv(pib_monthly,         "pib_monthly_vale_mucuri_2012_2023_v2.csv")
  readr::write_csv(check_benchmark,     "check_benchmark_vale_mucuri_2012_2023_v2.csv")
  readr::write_csv(teofilo_outliers,    "diagnostico_teofilo_outliers.csv")
  readr::write_csv(teofilo_seasonality, "diagnostico_teofilo_sazonalidade.csv")

  if (run_chow_lin) {
    readr::write_csv(chow_lin_monthly,  "pib_monthly_vale_mucuri_chow_lin.csv")
    readr::write_csv(chow_lin_check,    "check_benchmark_chow_lin.csv")
  }
}

# =========================================================
# 16) GRAFICOS RAPIDOS
# =========================================================

# 16.1 Ajuste anual in-sample
g_annual_fit <- panel_annual %>%
  ggplot(aes(x = pib_annual, y = pib_hat_annual_in)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(
    title = "Ajuste anual in-sample",
    x = "PIB anual oficial",
    y = "PIB anual ajustado pelo modelo"
  )

print(g_annual_fit)

# 16.2 Serie mensal benchmarked, Teofilo Otoni
g_teofilo_month <- pib_monthly %>%
  filter(norm_txt(name_muni) == "TEOFILO OTONI") %>%
  ggplot(aes(x = date, y = pib_month_est)) +
  geom_line() +
  labs(
    title = "PIB mensal benchmarked, Teofilo Otoni",
    x = NULL,
    y = "PIB mensal estimado"
  )

print(g_teofilo_month)

# 16.3 Sinal mensal limpo, Teofilo Otoni
g_teofilo_signal <- pib_monthly %>%
  filter(norm_txt(name_muni) == "TEOFILO OTONI") %>%
  ggplot(aes(x = date, y = signal_smooth)) +
  geom_line() +
  labs(
    title = "Sinal mensal limpo, Teofilo Otoni",
    x = NULL,
    y = "Signal smooth"
  )

print(g_teofilo_signal)

# 16.4 Cobertura media mensal, Teofilo Otoni
g_teofilo_cov <- pib_monthly %>%
  filter(norm_txt(name_muni) == "TEOFILO OTONI") %>%
  ggplot(aes(x = date, y = cf_cvg_mean)) +
  geom_line() +
  labs(
    title = "Cobertura media mensal, Teofilo Otoni",
    x = NULL,
    y = "cf_cvg_mean"
  )

print(g_teofilo_cov)

# =========================================================
# 17) RESUMOS FINAIS NO CONSOLE
# =========================================================
cat("\n=====================================================\n")
cat("RESUMO FINAL\n")
cat("=====================================================\n")
cat("Municipios:", length(unique(pib_monthly$code_muni)), "\n")
cat("Periodo mensal:", as.character(min(pib_monthly$date)), "a", as.character(max(pib_monthly$date)), "\n")
cat("Anos anuais oficiais:", min(panel_annual$year), "a", max(panel_annual$year), "\n")
cat("MAPE anual medio:", round(annual_metrics$mape, 3), "\n")
cat("Correlacao anual in-sample:", round(annual_metrics$cor, 4), "\n")
cat("Max erro de reconciliacao anual:", round(max(abs(check_benchmark$diff), na.rm = TRUE), 6), "\n")

if (run_chow_lin) {
  cat("Chow-Lin: executado.\n")
} else {
  cat("Chow-Lin: NAO executado. Altere run_chow_lin <- TRUE para ativar.\n")
}
cat("=====================================================\n")

# =========================================================
# FIM
# =========================================================
