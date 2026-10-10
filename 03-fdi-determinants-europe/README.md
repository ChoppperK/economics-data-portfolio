# Determinants of Foreign Direct Investment in Europe

*Economics of International Trade and Investment, Mae Fah Luang University · April 2025 · Python, Excel*

## Question

What economic factors attract foreign direct investment (FDI) into European countries?

## Data

A panel of **5 countries** (Denmark, Finland, France, Germany, Greece) over **2000–2023**, giving **120 observations**. Sources: UNCTAD and World Bank. I tested 10 explanatory variables, including imports, personal remittances, net migration, inflation, unemployment, GDP per capita, industry growth, and youth labor-force participation.

Files: [`data/fdi_stock_panel.csv`](data/fdi_stock_panel.csv), [`data/fdi_flow_panel.csv`](data/fdi_flow_panel.csv)

## Method

Two OLS regressions: one explaining the **stock** of FDI (total investment built up) and one explaining yearly **flows** (new investment). Full output: [`regression_results.txt`](regression_results.txt)

## Results

| | FDI stock model | FDI flow model |
|---|---|---|
| R² | **0.975** | 0.141 |
| Significant at 5% | Imports (+), remittances (+), net migration (−), inflation (−), unemployment (−) | GDP per capita (−) |

- **FDI stock is well explained** by the economy's size and openness: countries that import more and receive more remittances hold more foreign investment, while higher inflation and unemployment discourage it.
- **Yearly FDI flows are hard to predict** (R² = 0.14). New investment jumps around from year to year, often because of one-off deals.
- **Limitation:** with 5 countries the sample is small, and the very high condition number suggests multicollinearity among the variables.
