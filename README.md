# viirs-ntl-municipal-mg

Pipeline de prova de conceito para estimação de **proxies mensais de atividade econômica** para municípios brasileiros a partir de dados de **luzes noturnas VIIRS/Suomi NPP**, ancoradas no PIB municipal anual oficial do IBGE.

Este repositório contém o script R desenvolvido e validado para os **23 municípios da mesorregião Vale do Mucuri (Minas Gerais)**, cobrindo o período de **2012 a 2023**.

---

## O que este repositório faz

O script implementa um pipeline completo que:

1. Baixa as geometrias municipais do IBGE via pacote `geobr`
2. Extrai estatísticas mensais do VIIRS DNB (radiância + cobertura livre de nuvens) do Google Earth Engine via `rgee`
3. Baixa o PIB municipal anual oficial do IBGE/SIDRA via pacote `sidrar`
4. Aplica filtro de qualidade, winsorização e suavização ao sinal de luminosidade
5. Estima um modelo anual de calibração: `log(PIB) ~ log(luz) + efeitos fixos municipais + tendência` via `fixest`
6. Constrói uma proxy mensal benchmarked aos totais anuais oficiais (desagregação proporcional)
7. Implementa desagregação temporal via **Chow-Lin** (pacote `tempdisagg`)
8. Valida os resultados com métricas RMSE, MAPE e checagens de reconciliação anual
9. Gera diagnósticos automáticos e visualizações

---

## Principais características metodológicas

- **Sinal NTL ponderado por qualidade:** a radiância mensal é penalizada pela cobertura livre de nuvens relativa à mediana histórica de cada município, reduzindo o ruído de meses com baixa cobertura do satélite
- **Winsorização:** valores extremos de radiância são truncados nos percentis 1% e 99% da distribuição histórica de cada município
- **Média móvel centrada de 3 meses** para suavização do sinal
- **Desagregação temporal Chow-Lin** com o sinal NTL como indicador de alta frequência, garantindo que as estimativas mensais somem exatamente ao PIB anual oficial
- **Fallback automático** para desagregação proporcional quando o Chow-Lin não converge
- **Checagens de reconciliação anual** para verificar a consistência das estimativas mensais com os benchmarks anuais oficiais

---

## Estrutura do repositório

```
viirs-ntl-municipal-mg/
│
├── viirs_pib_vale_mucuri.R   # Script principal (pipeline completo)
└── README.md                 # Este arquivo
```

Arquivos de saída (gerados pelo script, não versionados):

| Arquivo | Descrição |
|---------|-----------|
| `lights_monthly_vale_mucuri_2012_2023_clean.csv` | Painel mensal de NTL com indicadores de qualidade |
| `panel_annual_vale_mucuri_2012_2023.csv` | Painel anual NTL + PIB com previsões do modelo |
| `pib_monthly_vale_mucuri_2012_2023_v2.csv` | Proxy mensal de PIB (método proporcional) |
| `pib_monthly_vale_mucuri_chow_lin.csv` | Proxy mensal de PIB (método Chow-Lin) |
| `check_benchmark_vale_mucuri_2012_2023_v2.csv` | Diagnósticos de reconciliação anual |
| `check_benchmark_chow_lin.csv` | Diagnósticos de reconciliação Chow-Lin |
| `diagnostico_teofilo_outliers.csv` | Diagnóstico de outliers para Teófilo Otoni |
| `diagnostico_teofilo_sazonalidade.csv` | Diagnóstico de sazonalidade para Teófilo Otoni |

---

## Requisitos

### Pacotes R

```r
pkgs <- c(
  "tidyverse", "lubridate", "sf", "geobr", "rgee", "sidrar",
  "fixest", "slider", "janitor", "stringi", "tempdisagg",
  "geojsonio", "zoo"
)
install.packages(pkgs)
```

### Google Earth Engine

Este script requer uma conta no Google Earth Engine e um projeto GEE configurado.

1. Crie uma conta gratuita em [https://earthengine.google.com/](https://earthengine.google.com/)
2. Crie um projeto GEE em [https://console.cloud.google.com/](https://console.cloud.google.com/)
3. Instale e configure o pacote `rgee` seguindo a documentação oficial: [https://r-spatial.github.io/rgee/](https://r-spatial.github.io/rgee/)
4. Autentique com `ee_install()` e `ee_Initialize()`

### Ambiente Python

O pacote `rgee` requer um ambiente Python com a API do Earth Engine instalada. Siga o guia de configuração do `rgee` para o seu sistema operacional.

---

## Como usar

1. Clone ou baixe este repositório
2. Abra `viirs_pib_vale_mucuri.R` no RStudio
3. Edite a seção de **Parâmetros Gerais** (Bloco 0) no início do script:

```r
# Caminho para o Python do seu ambiente rgee
reticulate_python <- "caminho/para/seu/python/rgee"

# ID do seu projeto no Google Earth Engine
ee_project <- "seu-projeto-gee"
```

4. No Bloco 4, substitua o e-mail placeholder pelo e-mail associado à sua conta GEE:

```r
ensure_rgee_sessioninfo(email = "seu-email@gmail.com")
```

5. Execute o script. O pipeline completo — incluindo a desagregação Chow-Lin — leva aproximadamente 10 a 20 minutos, dependendo da conexão com a internet e do tempo de resposta do GEE.

---

## Fontes de dados

| Fonte | Descrição | Acesso |
|-------|-----------|--------|
| NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG | Compostos mensais VIIRS Day/Night Band | Google Earth Engine |
| IBGE/SIDRA — Tabela 5938 | PIB dos Municípios, 2012–2023 | API pública via `sidrar` |
| IBGE — Malha municipal | Geometrias municipais, edição 2020 | Público via `geobr` |

---

## Referências

Beyer, R. C. M., Hu, Y., & Yao, J. (2022). *Measuring Quarterly Economic Growth from Outer Space*. IMF Working Papers, WP/22/109.

Gibson, J., Olivia, S., Boe-Gibson, G., & Li, C. (2021). Which night lights data should we use in economics, and where? *Journal of Development Economics*, 149, 102602.

Henderson, J. V., Storeygard, A., & Weil, D. N. (2012). Measuring Economic Growth from Outer Space. *American Economic Review*, 102(2), 994–1028.

---

## Licença

Este projeto é distribuído sob a [Licença MIT](https://opensource.org/licenses/MIT). Você é livre para usar, adaptar e distribuir o código com atribuição.

---

## Observações

- A série mensal de PIB produzida por este script é uma **proxy modelada**, não uma estatística oficial. Ela é ancorada nos valores anuais oficiais do IBGE, mas deve ser interpretada com a incerteza inerente a qualquer método de estimação indireta.
- O pipeline foi desenvolvido e validado para os municípios do Vale do Mucuri. O desempenho pode variar em municípios com baixo grau de urbanização ou predominância de atividade agrícola, onde o sinal NTL é menos informativo.
- Este repositório representa uma **aplicação piloto**. A metodologia foi projetada para escalar para os 853 municípios de Minas Gerais com adaptações menores.
