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
*work in progress*
