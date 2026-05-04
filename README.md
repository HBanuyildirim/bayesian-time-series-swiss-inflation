# bayesian-time-series-swiss-inflation
Bayesian time series analysis of Swiss inflation and consumer behavior with pre- and post-2020 regime comparison

## Research Question

Is Switzerland really protected from inflation, or does inflation still affect consumer behavior?

We test whether the relationship between inflation and retail consumption changes after 2020.

## Model

We estimate a Bayesian autoregressive model:

Retail growth_t = β1 + β2 Retail growth_{t-1} + β3 Inflation_{t-1} + β4 Interest rate_{t-1} + β5 FX growth_{t-1} + ε_t

- Bayesian estimation with conjugate priors
- Posterior simulation via Monte Carlo
- Comparison across regimes (pre- and post-2020)

## Data

- Swiss CPI (inflation)
- Retail sales (consumer demand)
- SNB interest rate
- CHF/USD exchange rate

Source: FRED database  
Period: 2000–2024  
Frequency: Quarterly


## Key Findings

- Before 2020:
  Inflation has a negative and statistically meaningful effect on retail consumption.

- After 2020:
  The relationship becomes highly uncertain and unstable.

- Interpretation:
  The results do not show a clear structural break, but rather an increase in uncertainty.

- Robustness:
  OLS results confirm the same pattern.

## Robustness Checks

- OLS regression
- Correlation analysis
- Comparison of coefficients across regimes

Results are consistent with Bayesian findings.

## Limitations

- Small sample size after 2020
- High uncertainty in posterior distributions
- Results should be interpreted with caution

## Contribution

This project shows how macroeconomic relationships can become unstable during periods of economic stress.

It combines:
- Bayesian time series modeling
- Real-world macroeconomic data
- Regime comparison analysis
## Repository Structure

- /figures → all plots
- /outputs → posterior summaries
- main.R → full analysis code

## Author

Master’s student in Finance & Money  
University of Basel  

Focus: Data-driven financial analysis, time series modeling, and macroeconomics
