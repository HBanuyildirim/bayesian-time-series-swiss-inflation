#install.packages(c("tidyverse", "readxl", "lubridate", "zoo", "ggplot2", "coda", "patchwork"))
install.packages("coda", dependencies = TRUE)
install.packages("patchwork")
if (any(!installed)) {
  install.packages(packages[!installed])
}
install.packages("fredr", dependencies = TRUE)

library(tidyverse)
library(fredr)
library(lubridate)
library(zoo)
library(coda)
library(ggplot2)

set.seed(123)

# ------------------------------------------------------------
# 1. FRED API key
# ------------------------------------------------------------
# 1) Go to FRED website
# 2) Create a free account
# 3) Get your API key
# 4) Replace YOUR_API_KEY_HERE below

fredr_set_key("Your Fred KEY")

# ------------------------------------------------------------
# 2. Download data automatically from FRED
# ------------------------------------------------------------

start_date <- as.Date("2000-01-01")

retail_raw <- fredr(
  series_id = "SLRTTO02CHQ661N",
  observation_start = start_date
)

cpi_raw <- fredr(
  series_id = "CHECPIALLMINMEI",
  observation_start = start_date
)

rate_raw <- fredr(
  series_id = "IR3TIB01CHM156N",
  observation_start = start_date
)

fx_raw <- fredr(
  series_id = "DEXSZUS",
  observation_start = start_date
)

# ------------------------------------------------------------
# 3. Convert all data to quarterly frequency
# ------------------------------------------------------------

retail_q <- retail_raw %>%
  transmute(
    Date = as.yearqtr(date),
    RetailSales = value
  )

cpi_q <- cpi_raw %>%
  transmute(
    Date = as.yearqtr(date),
    CPI = value
  ) %>%
  group_by(Date) %>%
  summarise(CPI = mean(CPI, na.rm = TRUE), .groups = "drop")

rate_q <- rate_raw %>%
  transmute(
    Date = as.yearqtr(date),
    Rate = value
  ) %>%
  group_by(Date) %>%
  summarise(Rate = mean(Rate, na.rm = TRUE), .groups = "drop")

fx_q <- fx_raw %>%
  transmute(
    Date = as.yearqtr(date),
    FX = value
  ) %>%
  group_by(Date) %>%
  summarise(FX = mean(FX, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 4. Merge and prepare variables
# ------------------------------------------------------------

df <- retail_q %>%
  left_join(cpi_q, by = "Date") %>%
  left_join(rate_q, by = "Date") %>%
  left_join(fx_q, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    Date_real = as.Date(Date),
    
    # Retail sales growth
    # Similar to GDP growth in Assignment 1.3:
    # Δy_t = 100 * (log(Y_t) - log(Y_{t-1}))
    retail_growth = 100 * (log(RetailSales) - log(lag(RetailSales, 1))),
    
    # Inflation
    inflation = 100 * (log(CPI) - log(lag(CPI, 1))),
    
    # Exchange rate growth
    fx_growth = 100 * (log(FX) - log(lag(FX, 1))),
    
    # Lagged explanatory variables
    retail_growth_lag = lag(retail_growth, 1),
    inflation_lag = lag(inflation, 1),
    rate_lag = lag(Rate, 1),
    fx_growth_lag = lag(fx_growth, 1)
  ) %>%
  drop_na()

# ------------------------------------------------------------
# 5. Plot raw transformed data
# ------------------------------------------------------------

p_data <- df %>%
  select(Date_real, retail_growth, inflation, Rate, fx_growth) %>%
  pivot_longer(
    cols = c(retail_growth, inflation, Rate, fx_growth),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  ggplot(aes(x = Date_real, y = Value)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = "Swiss Macroeconomic Indicators",
    subtitle = "Retail sales growth, inflation, interest rate, and CHF/USD exchange rate growth",
    x = "",
    y = ""
  )

print(p_data)

# ------------------------------------------------------------
# 6. Bayesian regression function
# Based on Assignment 1.3 logic:
#
# y_t = beta_1 + beta_2 y_{t-1} + beta_3 x_{t-1} + error_t
#
# Here:
# retail_growth_t =
# beta_1
# + beta_2 retail_growth_{t-1}
# + beta_3 inflation_{t-1}
# + beta_4 rate_{t-1}
# + beta_5 fx_growth_{t-1}
# + error_t
# ------------------------------------------------------------

bayesian_regression <- function(data_input, M = 5000) {
  
  y <- as.matrix(data_input$retail_growth)
  
  X <- data_input %>%
    select(
      retail_growth_lag,
      inflation_lag,
      rate_lag,
      fx_growth_lag
    ) %>%
    as.matrix()
  
  X <- cbind(1, X)
  
  T_obs <- nrow(X)
  k <- ncol(X)
  
  # Prior:
  # beta | sigma^2 ~ N(b0, sigma^2 B0)
  # sigma^2 ~ IG(g0, G0)
  #
  # Same idea as Assignment 1.3:
  # b0 = 0
  # B0 = T * (X'X)^(-1)
  
  b0 <- matrix(0, k, 1)
  B0 <- T_obs * solve(t(X) %*% X)
  B0_inv <- solve(B0)
  
  g0 <- 2
  G0 <- 1
  
  # Posterior moments
  B1 <- solve(t(X) %*% X + B0_inv)
  b1 <- B1 %*% (t(X) %*% y + B0_inv %*% b0)
  
  g1 <- g0 + T_obs / 2
  
  G1 <- G0 + 0.5 * (
    t(y) %*% y +
      t(b0) %*% B0_inv %*% b0 -
      t(b1) %*% solve(B1) %*% b1
  )
  
  G1 <- as.numeric(G1)
  
  # Draw sigma^2 from inverse-gamma
  sigma2_draws <- 1 / rgamma(M, shape = g1, rate = G1)
  
  # Draw beta conditional on sigma^2
  beta_draws <- matrix(NA, nrow = M, ncol = k)
  
  for (m in 1:M) {
    beta_draws[m, ] <- as.numeric(
      b1 + t(chol(sigma2_draws[m] * B1)) %*% rnorm(k)
    )
  }
  
  colnames(beta_draws) <- c(
    "Intercept",
    "Lagged retail growth",
    "Lagged inflation",
    "Lagged interest rate",
    "Lagged FX growth"
  )
  
  return(list(
    beta_draws = beta_draws,
    sigma2_draws = sigma2_draws,
    posterior_mean = b1,
    posterior_variance = B1
  ))
}

# ------------------------------------------------------------
# 7. Estimate full sample, pre-2020, and post-2020 models
# ------------------------------------------------------------

full_model <- bayesian_regression(df, M = 5000)

df_pre2020 <- df %>%
  filter(Date_real < as.Date("2020-01-01"))

df_post2020 <- df %>%
  filter(Date_real >= as.Date("2020-01-01"))

pre_model <- bayesian_regression(df_pre2020, M = 5000)
post_model <- bayesian_regression(df_post2020, M = 5000)

# ------------------------------------------------------------
# 8. Posterior summary function
# ------------------------------------------------------------

posterior_summary <- function(model) {
  
  beta <- model$beta_draws
  sigma <- model$sigma2_draws
  
  beta_summary <- data.frame(
    Parameter = colnames(beta),
    Mean = apply(beta, 2, mean),
    Variance = apply(beta, 2, var),
    HPDI_Lower = apply(beta, 2, function(x) {
      HPDinterval(as.mcmc(x), prob = 0.95)[1]
    }),
    HPDI_Upper = apply(beta, 2, function(x) {
      HPDinterval(as.mcmc(x), prob = 0.95)[2]
    })
  )
  
  sigma_summary <- data.frame(
    Parameter = "Sigma2",
    Mean = mean(sigma),
    Variance = var(sigma),
    HPDI_Lower = HPDinterval(as.mcmc(sigma), prob = 0.95)[1],
    HPDI_Upper = HPDinterval(as.mcmc(sigma), prob = 0.95)[2]
  )
  
  bind_rows(beta_summary, sigma_summary)
}

full_summary <- posterior_summary(full_model)
pre_summary <- posterior_summary(pre_model)
post_summary <- posterior_summary(post_model)

print("Full sample posterior summary")
print(full_summary)

print("Pre-2020 posterior summary")
print(pre_summary)

print("Post-2020 posterior summary")
print(post_summary)

# ------------------------------------------------------------
# 9. Plot posterior densities
# ------------------------------------------------------------

plot_posterior <- function(model, title_text) {
  
  beta_df <- as.data.frame(model$beta_draws) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Draw"
    )
  
  ggplot(beta_df, aes(x = Draw)) +
    geom_density(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ Parameter, scales = "free") +
    theme_minimal() +
    labs(
      title = title_text,
      x = "Posterior draw",
      y = "Density"
    )
}

p_full <- plot_posterior(full_model, "Posterior Distributions - Full Sample")
p_pre <- plot_posterior(pre_model, "Posterior Distributions - Pre-2020")
p_post <- plot_posterior(post_model, "Posterior Distributions - Post-2020")

print(p_full)
print(p_pre)
print(p_post)

# ------------------------------------------------------------
# 10. Compare inflation effect before and after 2020
# ------------------------------------------------------------

inflation_compare <- data.frame(
  Period = c(
    rep("Pre-2020", nrow(pre_model$beta_draws)),
    rep("Post-2020", nrow(post_model$beta_draws))
  ),
  Inflation_Effect = c(
    pre_model$beta_draws[, "Lagged inflation"],
    post_model$beta_draws[, "Lagged inflation"]
  )
)

p_inflation_compare <- ggplot(
  inflation_compare,
  aes(x = Inflation_Effect, fill = Period)
) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Does the Inflation Effect Change After 2020?",
    subtitle = "Posterior draws of the lagged inflation coefficient",
    x = "Inflation coefficient",
    y = "Density"
  )

print(p_inflation_compare)

# ------------------------------------------------------------
# 11. Compare all coefficients pre vs post 2020
# ------------------------------------------------------------

coef_compare <- bind_rows(
  as.data.frame(pre_model$beta_draws) %>%
    mutate(Period = "Pre-2020"),
  
  as.data.frame(post_model$beta_draws) %>%
    mutate(Period = "Post-2020")
) %>%
  pivot_longer(
    cols = -Period,
    names_to = "Parameter",
    values_to = "Draw"
  )

p_coef_compare <- ggplot(coef_compare, aes(x = Draw, fill = Period)) +
  geom_density(alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(
    title = "Pre-2020 vs Post-2020 Posterior Comparison",
    x = "Posterior draw",
    y = "Density"
  )

print(p_coef_compare)

# ------------------------------------------------------------
# 12. Save outputs for GitHub
# ------------------------------------------------------------

if (!dir.exists("figures")) {
  dir.create("figures")
}

if (!dir.exists("outputs")) {
  dir.create("outputs")
}

ggsave("figures/01_swiss_macro_indicators.png", p_data, width = 10, height = 7)
ggsave("figures/02_posterior_full_sample.png", p_full, width = 10, height = 6)
ggsave("figures/03_posterior_pre2020.png", p_pre, width = 10, height = 6)
ggsave("figures/04_posterior_post2020.png", p_post, width = 10, height = 6)
ggsave("figures/05_inflation_effect_pre_post_2020.png", p_inflation_compare, width = 8, height = 5)
ggsave("figures/06_all_coefficients_pre_post_2020.png", p_coef_compare, width = 10, height = 6)

write_csv(full_summary, "outputs/posterior_summary_full.csv")
write_csv(pre_summary, "outputs/posterior_summary_pre2020.csv")
write_csv(post_summary, "outputs/posterior_summary_post2020.csv")

# ------------------------------------------------------------
# 13. Simple interpretation helper
# ------------------------------------------------------------

interpret_parameter <- function(summary_table, parameter_name, period_name) {
  
  row <- summary_table %>%
    filter(Parameter == parameter_name)
  
  if (nrow(row) == 0) {
    print(paste("Parameter not found:", parameter_name))
    return(NULL)
  }
  
  if (row$HPDI_Lower > 0 | row$HPDI_Upper < 0) {
    print(paste(
      period_name, "-",
      parameter_name,
      "is relevant: 95% HPDI does not include zero."
    ))
  } else {
    print(paste(
      period_name, "-",
      parameter_name,
      "is not clearly different from zero: 95% HPDI includes zero."
    ))
  }
}

interpret_parameter(pre_summary, "Lagged inflation", "Pre-2020")
interpret_parameter(post_summary, "Lagged inflation", "Post-2020")

# ============================================================
# End of script
# ============================================================

lm_full <- lm(retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag, data = df)
summary(lm_full)

lm_pre <- lm(retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag, data = df_pre2020)
summary(lm_pre)

lm_post <- lm(retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag, data = df_post2020)
summary(lm_post)


cor(df_pre2020$retail_growth, df_pre2020$inflation_lag, use = "complete.obs")
cor(df_post2020$retail_growth, df_post2020$inflation_lag, use = "complete.obs")

retail_growth = 100 * (log(RetailSales) - log(lag(RetailSales, 1)))
inflation = 100 * (log(CPI) - log(lag(CPI, 1)))

retail_growth = 100 * (log(RetailSales) - log(lag(RetailSales, 4)))
inflation = 100 * (log(CPI) - log(lag(CPI, 4)))
fx_growth = 100 * (log(FX) - log(lag(FX, 4)))



# ============================================================
# Robustness Check: OLS and Simple Correlation
# ============================================================

# 1. Check column names
colnames(df)

# 2. OLS models
lm_full <- lm(
  retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag,
  data = df
)

lm_pre <- lm(
  retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag,
  data = df_pre2020
)

lm_post <- lm(
  retail_growth ~ retail_growth_lag + inflation_lag + rate_lag + fx_growth_lag,
  data = df_post2020
)

summary(lm_full)
summary(lm_pre)
summary(lm_post)

# 3. Simple correlation check
cor_pre <- cor(
  df_pre2020$retail_growth,
  df_pre2020$inflation_lag,
  use = "complete.obs"
)

cor_post <- cor(
  df_post2020$retail_growth,
  df_post2020$inflation_lag,
  use = "complete.obs"
)

print(paste("Pre-2020 correlation between retail growth and lagged inflation:", round(cor_pre, 3)))
print(paste("Post-2020 correlation between retail growth and lagged inflation:", round(cor_post, 3)))

# 4. Compare OLS coefficients directly
ols_compare <- data.frame(
  Parameter = names(coef(lm_pre)),
  Pre_2020 = coef(lm_pre),
  Post_2020 = coef(lm_post)
)

print(ols_compare)