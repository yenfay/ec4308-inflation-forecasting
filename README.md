# EC4308 Project: Inflation Forecasting
By: 

## Motivation and Research Question
Accurate inflation forecasting is vital for effective and rational economic decision-making. Central banks rely on inflation projections to set expectations and interest rates, while firms and households depend on these projections for long-term contract, wage, and investment decision-making. Traditional linear models often struggle to maintain predictive accuracy during volatile economic regimes, whereas modern machine learning (ML) estimators offer greater flexibility in handling high-dimensional macroeconomic panels.

**Research question:** With information on high-dimensional macroeconomic indicators, how effective are ML approaches in forecasting Canadian inflation, and how does their predictive performance vary across forecast horizons?

## Data Sources and Preprocessing
We use the stationary and balanced monthly panel (balanced_can_md) of the LCDMA dataset (Fortin-Gagnon et al., 2022), which begins in 1981M01 and contains 410 national and provincial Canadian macroeconomic and financial indicators. It is already transformed to stationarity, with missing observations imputed and relevant series seasonally adjusted. As CPI is provided without seasonal adjustment, we construct the target from Statistics Canada’s seasonally adjusted CPI. As secondary analysis, we add a few US macroeconomic indicators from FRED-MD (McCracken and Ng, 2016) to account for global economic influences; these are already seasonally adjusted and are transformed in the same way as their Canadian counterparts.

## Methodology and Forecast Evaluation
Subject to computational constraints, we compare a diverse set of models, from penalised regressions to ensemble methods, to capture different properties of the macroeconomic signals:
1. Traditional benchmarks: AR(p) and Random Walk.
2. Penalized regressions: LASSO, Post-LASSO, and Ridge.
3. Dimension reduction: Principal Component Regression and Partial Least Squares.
4. Non-linear / ensemble: Random Forest and XGBoost.
5. Hybrid: Factor-Augmented Sparse Linear Regression Model (Fan et al., 2024).

Forecasts are generated with a rolling window, allowing models adapt to structural breaks and changing regimes. Given temporal dependence, model parameters are tuned by rolling-origin cross-validation strictly within each training window. Forecast accuracy is evaluated under MSE loss, and the Diebold–Mariano (DM) test assesses whether squared-error loss differentials between models are significant; the rolling-window scheme resolves the nested-model issue.

## Repository Structure

```
ec4308-inflation-forecasting/
├── data/
│   ├── raw/                  LCDMA balanced_can_md, FRED-MD (2026-rev-08-MD), StatCan CPI tables 18-10-0004 / 18-10-0006
│   └── processed/            {main,extended}_data_{self,official}.rds (cleaning.Rmd output)
├── cleaning.Rmd              raw → processed panels (CPI target construction, FRED-MD US extension)
├── eda.Rmd                   exploratory analysis of processed panels (no files written)
├── benchmarks.R              Random Walk and AR(p) (BIC) forecasts
├── models/
│   ├── lasso.Rmd             LASSO / Post-LASSO (sources func-lasso.R); one panel per run, set in config
│   ├── func-lasso.R          rolling-window LASSO / Post-LASSO helpers
│   ├── ridge_farm.Rmd        Ridge (+ 4 PCs) and FARM, incl. US-extension and US-key-indicator runs
│   └── rf-xgb-test.Rmd       Random Forest, Bagging and XGBoost
├── evaluate.Rmd              forecast combinations, MSE / relative MSE, DM tests, robustness, plots
└── results/
    ├── forecasts/            benchmark and RF/Bag/XGB forecasts, AR lag choices, RF/XGB tuning paths
    ├── lasso/                LASSO forecasts (.rds + .csv), tuning paths, fitted result objects
    ├── farm-ridge/           Ridge / FARM forecasts (main, _us, _uskey)
    └── evaluation/           metrics, DM tables, combination forecasts, figures
```

### Pipeline

```
 data/raw ──► cleaning.Rmd ──► data/processed/*.rds
                                   │
        ┌──────────────┬───────────┼────────────────┬──────────────────┐
        ▼              ▼           ▼                ▼                  ▼
     eda.Rmd     benchmarks.R  models/lasso.Rmd  models/ridge_      models/rf-xgb-
  (diagnostics)    (RW, AR)    (+ func-lasso.R)  farm.Rmd           test.Rmd
                       │      main + extended    main + extended    (RF, Bag, XGB)
                    main only        │                │              main only
                       ▼             ▼                ▼                  ▼
              results/forecasts  results/lasso  results/farm-ridge  results/forecasts
                       └─────────────┴───────┬────────┴──────────────────┘
                                             ▼
                                        evaluate.Rmd
                                             ▼
                                    results/evaluation/
              (combinations, MSE by horizon, DM tests, official-vs-self, US extension)
```

**Inputs by script**

| Script | Panels read | Horizons | Outputs |
|---|---|---|---|
| `benchmarks.R` | `main_data_{official,self}` | 1, 3, 6, 12 | `results/forecasts/forecasts_benchmarks_<sample>.csv`, `ar_selected_lags_<sample>.csv` |
| `models/lasso.Rmd` | one of the four panels per run (`DATASET`, optional `USE_US_KEY` basket) | 1, 3, 6, 12 | `results/lasso/forecasts_<run>.rds`, `results_lasso_<run>.rds`; then `forecasts_lasso_<sample>.csv` and `tuning_path_lasso_<sample>.csv` for every run whose `.rds` exists |
| `models/ridge_farm.Rmd` | `main_data_*` (main run), `extended_data_*` (all US series, and the US-key basket) | 1–12 | `results/farm-ridge/forecasts_ridge_farm_<sample>[_us\|_uskey].csv`, tuning-path CSVs; working `.rds` files in `results/` |
| `models/rf-xgb-test.Rmd` | `main_data_{official,self}` | 1, 3, 6, 12 | `results/forecasts/forecasts_rf_xgb_<sample>.csv`, `tuning_path_rf_xgb_<sample>.csv` |

**Running order.** Run `cleaning.Rmd` first, then the model scripts in any order, then `evaluate.Rmd`. `benchmarks.R` and `evaluate.Rmd` use paths relative to the repository root, so run them with the working directory at the root; the `models/` scripts resolve paths with `here::here()`. `lasso.Rmd` is run once per panel: set `DATASET` / `USE_US_KEY` in the config chunk (e.g. `main_data_official`, `main_data_self`, `extended_data_official`, `extended_data_self`), and the export chunks then write CSVs for every run that has been saved. `ridge_farm.Rmd` loads saved `.rds` results by default (`MODELS_TO_RUN <- character(0)`); set `MODELS_TO_RUN <- c("ridge", "farm")` to re-estimate.

**What `evaluate.Rmd` consumes.**
- Main comparison (models RW, AR_BIC, LASSO, Post-LASSO, Ridge, FARM, RF, Bagging, XGB): `forecasts_benchmarks_<sample>.csv`, `forecasts_ridge_farm_<sample>.csv`, `forecasts_rf_xgb_<sample>.csv` and the LASSO `.rds` `results/lasso/forecasts_main_data_<sample>.rds`. The official sample is required; the self sample is optional and skipped if any file is missing.
- All models are evaluated on a common window, target dates 2012-02 to 2026-07 (174 targets), at h = 1, 3, 6, 12, for both CPI targets. Self-sample forecasts for target dates before 2012-02 are dropped, as are horizons other than 1, 3, 6 and 12 (e.g. Ridge/FARM's 1–12).
- US extension (optional, used if present): `forecasts_ridge_farm_official_us.csv`, `forecasts_ridge_farm_official_uskey.csv`, `forecasts_lasso_official_us.csv`, `forecasts_lasso_self_us.csv`. The `*_self_us`/`*_self_uskey` Ridge/FARM files and the LASSO `*_uskey` files are produced but not currently evaluated.
