# =============================================================================
# SUR LECTURE: WHEN CAN SEPARATE EQUATIONS HELP ESTIMATE EACH OTHER?
# =============================================================================
#
# Purpose
# -------
# Reproduce the live econometrics lab for the General Electric (GE) and
# Westinghouse investment equations. The script deliberately calculates SUR
# from the GLS matrix formula so students can see how the estimator works.
#
# Run this script from its folder. The data are entered below, so no external
# data file or internet connection is required.
# =============================================================================

# ---- Packages ---------------------------------------------------------------

library(tidyverse)
library(broom)

# Save figures next to this script, even when it is launched from elsewhere.
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_directory <- if (length(script_argument) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_argument)))
} else {
  getwd()
}

figure_directory <- file.path(script_directory, "figures")
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

# ---- "Baby Judge" Textbook data, Table 17.1---------------------------------

# i = investment, v = market value, k = existing capital stock.
# The textbook labels the observations 1,...,20; no calendar-year mapping is
# required for the exercise.
investment <- tribble(
  ~year, ~i1,  ~v1,    ~k1,  ~i2,  ~v2,   ~k2,
      1, 33.1, 1170.6,  97.8, 12.93, 191.5,   1.8,
      2, 45.0, 2015.8, 104.4, 25.90, 516.0,   0.8,
      3, 77.2, 2803.3, 118.0, 35.05, 729.0,   7.4,
      4, 44.6, 2039.7, 156.2, 22.89, 560.4,  18.1,
      5, 48.1, 2256.2, 172.6, 18.84, 519.9,  23.5,
      6, 74.4, 2132.2, 186.6, 28.57, 628.5,  26.5,
      7,113.0, 1834.1, 220.9, 48.51, 537.1,  36.2,
      8, 91.9, 1588.0, 287.8, 43.34, 561.2,  60.8,
      9, 61.3, 1749.4, 319.9, 37.02, 617.2,  84.4,
     10, 56.8, 1687.2, 321.3, 37.81, 626.7,  91.2,
     11, 93.6, 2007.7, 319.6, 39.27, 737.2,  92.4,
     12,159.9, 2208.3, 346.0, 53.46, 760.5,  86.0,
     13,147.2, 1656.7, 456.4, 55.56, 581.4, 111.1,
     14,146.3, 1604.4, 543.4, 49.56, 662.3, 130.6,
     15, 98.3, 1431.8, 618.3, 32.04, 583.8, 141.8,
     16, 93.5, 1610.5, 647.4, 32.24, 635.2, 136.7,
     17,135.2, 1819.4, 671.3, 54.38, 723.8, 129.7,
     18,157.3, 2079.7, 726.1, 71.78, 864.1, 145.5,
     19,179.5, 2371.6, 800.3, 90.08,1193.5, 174.8,
     20,189.6, 2759.9, 888.9, 68.60,1188.9, 213.5
)

investment

# =============================================================================
# LAB STEP 1: ESTABLISH THE SEPARATE-EQUATION OLS BENCHMARK
# =============================================================================

# PREDICTION STOP:
# Before running the models, predict the signs of the market-value and capital
# coefficients. Which coefficients do you expect to be estimated most precisely?

ge_ols <- lm(i1 ~ v1 + k1, data = investment)
west_ols <- lm(i2 ~ v2 + k2, data = investment)

ge_ols_table <- tidy(ge_ols, conf.int = TRUE) %>%
  mutate(equation = "General Electric", .before = 1)

west_ols_table <- tidy(west_ols, conf.int = TRUE) %>%
  mutate(equation = "Westinghouse", .before = 1)

ols_table <- bind_rows(ge_ols_table, west_ols_table)

print(ols_table)

# Key benchmark estimates (rounded):
# GE:          intercept = -9.9563, v1 = 0.02655, k1 = 0.15169
# Westinghouse: intercept = -0.5094, v2 = 0.05289, k2 = 0.09241

# =============================================================================
# LAB STEP 2: DIAGNOSE WHETHER THE EQUATIONS SHARE UNEXPLAINED SHOCKS
# =============================================================================

# PREDICTION STOP:
# If both firms experience common macroeconomic and credit-market shocks, what
# sign should the contemporaneous residual correlation have?

residual_data <- investment %>%
  transmute(
    year,
    ge_residual = residuals(ge_ols),
    west_residual = residuals(west_ols)
  )

rho_hat <- cor(
  residual_data$ge_residual,
  residual_data$west_residual
)

# For M = 2 equations, the Breusch-Pagan LM statistic is T*rho_hat^2 and
# has an asymptotic chi-squared distribution with one degree of freedom.
T_obs <- nrow(investment)
lm_statistic <- T_obs * rho_hat^2
lm_p_value <- pchisq(lm_statistic, df = 1, lower.tail = FALSE)

diagnostic_results <- tibble(
  residual_correlation = rho_hat,
  LM_statistic = lm_statistic,
  degrees_of_freedom = 1,
  p_value = lm_p_value
)

print(diagnostic_results)

# Verified values:
# rho_hat = 0.728965; LM = 10.62780; p = 0.001114.

residual_long <- residual_data %>%
  pivot_longer(
    cols = ends_with("residual"),
    names_to = "firm",
    values_to = "residual"
  ) %>%
  mutate(
    firm = recode(
      firm,
      ge_residual = "General Electric",
      west_residual = "Westinghouse"
    )
  )

residual_series_plot <- ggplot(
  residual_long,
  aes(x = year, y = residual, colour = firm)
) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#2563A6", "#D97706")) +
  labs(
    title = "OLS residuals move together in many periods",
    subtitle = "Contemporaneous residual correlation = 0.729",
    x = "Observation",
    y = "OLS residual",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

residual_scatter_plot <- ggplot(
  residual_data,
  aes(x = ge_residual, y = west_residual)
) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 2.6, colour = "#17324D") +
  geom_smooth(method = "lm", se = FALSE, colour = "#D97706") +
  labs(
    title = "The equation errors are statistically related",
    subtitle = "A common shock is one interpretation—not a causal link",
    x = "GE OLS residual",
    y = "Westinghouse OLS residual"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(figure_directory, "residual-series.png"),
  residual_series_plot,
  width = 9,
  height = 5.2,
  dpi = 180
)

ggsave(
  file.path(figure_directory, "residual-scatter.png"),
  residual_scatter_plot,
  width = 7.2,
  height = 5.2,
  dpi = 180
)

# =============================================================================
# LAB STEP 3: DELIBERATE FAILURE -- STACK THE EQUATIONS AND APPLY OLS
# =============================================================================

# PLAUSIBLE BUT WRONG CLAIM:
# "Once the equations are stacked, OLS uses the cross-equation correlation."
#
# Because the system regressor matrix is block diagonal, system OLS is exactly
# the same estimator as applying OLS separately to each equation.

X1 <- model.matrix(ge_ols)
X2 <- model.matrix(west_ols)
y1 <- investment$i1
y2 <- investment$i2

zero_12 <- matrix(0, nrow = T_obs, ncol = ncol(X2))
zero_21 <- matrix(0, nrow = T_obs, ncol = ncol(X1))

X_system <- rbind(
  cbind(X1, zero_12),
  cbind(zero_21, X2)
)

y_system <- c(y1, y2)

beta_system_ols <- solve(
  crossprod(X_system),
  crossprod(X_system, y_system)
)

beta_separate_ols <- c(coef(ge_ols), coef(west_ols))
maximum_ols_difference <- max(abs(beta_system_ols - beta_separate_ols))

print(
  tibble(
    coefficient = c(
      "GE intercept", "GE market value", "GE capital",
      "Westinghouse intercept", "Westinghouse market value",
      "Westinghouse capital"
    ),
    separate_OLS = beta_separate_ols,
    system_OLS = as.numeric(beta_system_ols),
    difference = as.numeric(beta_system_ols) - beta_separate_ols
  )
)

cat("Maximum absolute OLS difference:", maximum_ols_difference, "\n")
stopifnot(maximum_ols_difference < 1e-8)

# =============================================================================
# LAB STEP 4: ESTIMATE FEASIBLE SUR FROM THE GLS MATRIX FORMULA
# =============================================================================

# Step 4a: Estimate the 2 x 2 contemporaneous covariance matrix from the two
# sets of OLS residuals. The equations each have three coefficients, so T - 3
# is the common degrees-of-freedom divisor used here.

E_hat <- cbind(
  GE = residuals(ge_ols),
  Westinghouse = residuals(west_ols)
)

K_average <- mean(c(ncol(X1), ncol(X2)))
Sigma_hat <- crossprod(E_hat) / (T_obs - K_average)

# Step 4b: Expand the 2 x 2 covariance matrix to the 40 x 40 system covariance
# matrix. The identity matrix represents no correlation across different years.
Omega_hat <- kronecker(Sigma_hat, diag(T_obs))
Omega_hat_inverse <- solve(Omega_hat)

# Step 4c: Apply the familiar GLS estimator to the stacked system.
sur_covariance <- solve(
  t(X_system) %*% Omega_hat_inverse %*% X_system
)

beta_sur <- sur_covariance %*%
  t(X_system) %*%
  Omega_hat_inverse %*%
  y_system

sur_standard_error <- sqrt(diag(sur_covariance))
sur_z_value <- as.numeric(beta_sur) / sur_standard_error
sur_p_value <- 2 * pnorm(abs(sur_z_value), lower.tail = FALSE)

sur_table <- tibble(
  equation = rep(c("General Electric", "Westinghouse"), each = 3),
  term = rep(c("Intercept", "Market value", "Capital"), times = 2),
  estimate = as.numeric(beta_sur),
  std_error = sur_standard_error,
  statistic = sur_z_value,
  p_value = sur_p_value
)

print(Sigma_hat)
print(sur_table)

# Verified SUR estimates:
# GE:          intercept = -27.7193, v1 = 0.03831, k1 = 0.13904
# Westinghouse: intercept = -1.2520, v2 = 0.05763, k2 = 0.06398

# Compare OLS and SUR slopes. Multiplying estimates and standard errors by 100
# makes the interpretation a 100-unit increase in market value or capital.
ols_slopes <- ols_table %>%
  filter(term != "(Intercept)") %>%
  transmute(
    equation,
    term = if_else(str_detect(term, "^v"), "Market value", "Capital"),
    estimator = "Separate OLS",
    estimate_100 = 100 * estimate,
    lower_100 = 100 * (estimate - 1.96 * std.error),
    upper_100 = 100 * (estimate + 1.96 * std.error)
  )

sur_slopes <- sur_table %>%
  filter(term != "Intercept") %>%
  transmute(
    equation,
    term,
    estimator = "SUR",
    estimate_100 = 100 * estimate,
    lower_100 = 100 * (estimate - 1.96 * std_error),
    upper_100 = 100 * (estimate + 1.96 * std_error)
  )

coefficient_comparison <- bind_rows(ols_slopes, sur_slopes) %>%
  mutate(
    estimator = factor(estimator, levels = c("Separate OLS", "SUR")),
    term = factor(term, levels = c("Market value", "Capital"))
  )

coefficient_plot <- ggplot(
  coefficient_comparison,
  aes(
    x = estimate_100,
    y = term,
    xmin = lower_100,
    xmax = upper_100,
    colour = estimator
  )
) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_errorbar(
    position = position_dodge(width = 0.45),
    height = 0,
    linewidth = 0.8
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.8) +
  facet_wrap(~equation, scales = "free_x") +
  scale_colour_manual(values = c("#6B7280", "#2563A6")) +
  labs(
    title = "SUR changes both estimates and estimated precision",
    subtitle = "Effects shown for a 100-unit increase; bars are 95% intervals",
    x = "Estimated change in investment",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(
  file.path(figure_directory, "coefficient-comparison.png"),
  coefficient_plot,
  width = 9,
  height = 5.3,
  dpi = 180
)

# =============================================================================
# LAB STEP 5: STRESS-TEST THE CLAIM THAT CORRELATION IS ENOUGH
# =============================================================================

# This is a mathematical demonstration, not a proposed economic model. Give
# both equations the same X matrix. SUR must then reproduce separate OLS even
# if their residuals are correlated.

X2_same <- X1
beta2_same_ols <- solve(crossprod(X2_same), crossprod(X2_same, y2))
e2_same <- y2 - X2_same %*% beta2_same_ols

E_same <- cbind(residuals(ge_ols), as.numeric(e2_same))
Sigma_same <- crossprod(E_same) / (T_obs - ncol(X1))

X_same_system <- rbind(
  cbind(X1, matrix(0, T_obs, ncol(X1))),
  cbind(matrix(0, T_obs, ncol(X1)), X1)
)

Omega_same <- kronecker(Sigma_same, diag(T_obs))
Omega_same_inverse <- solve(Omega_same)

beta_same_sur <- solve(
  t(X_same_system) %*% Omega_same_inverse %*% X_same_system,
  t(X_same_system) %*% Omega_same_inverse %*% y_system
)

beta_same_ols <- c(coef(ge_ols), beta2_same_ols)
identical_X_difference <- max(abs(beta_same_sur - beta_same_ols))

cat(
  "Maximum |SUR - OLS| when X1 = X2:",
  identical_X_difference,
  "\n"
)

stopifnot(identical_X_difference < 1e-8)

# ---- Exact efficiency map for a simplified two-equation model ---------------

# In a system with one standardized regressor per equation, equal error
# variances, error correlation rho_e, and regressor correlation r_x:
#
# Var(SUR) / Var(OLS) = (1 - rho_e^2) / (1 - rho_e^2 * r_x^2).
#
# A ratio of 1 means no gain. A smaller ratio means greater SUR precision.

efficiency_grid <- crossing(
  error_correlation = seq(0, 0.95, by = 0.01),
  regressor_overlap = seq(0, 1, by = 0.01)
) %>%
  mutate(
    relative_variance =
      (1 - error_correlation^2) /
      (1 - error_correlation^2 * regressor_overlap^2),
    precision_gain_percent = 100 * (1 - relative_variance)
  )

efficiency_plot <- ggplot(
  efficiency_grid,
  aes(
    x = regressor_overlap,
    y = error_correlation,
    fill = precision_gain_percent
  )
) +
  geom_raster() +
  geom_contour(
    aes(z = precision_gain_percent),
    breaks = c(10, 25, 50, 75),
    colour = "white",
    linewidth = 0.35
  ) +
  scale_fill_viridis_c(option = "C", limits = c(0, 90)) +
  labs(
    title = "SUR gains need correlated errors and different regressors",
    subtitle = "Exact result for a standardized two-equation, one-regressor model",
    x = expression("Regressor overlap  " * r[x]),
    y = expression("Error correlation  " * rho[e]),
    fill = "Variance\nreduction (%)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(figure_directory, "efficiency-map.png"),
  efficiency_plot,
  width = 8.4,
  height = 5.8,
  dpi = 180
)

# =============================================================================
# FINAL COMPARISON TABLE
# =============================================================================

ols_for_comparison <- ols_table %>%
  transmute(
    equation,
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      str_detect(term, "^v") ~ "Market value",
      TRUE ~ "Capital"
    ),
    OLS_estimate = estimate,
    OLS_SE = std.error
  )

final_comparison <- ols_for_comparison %>%
  left_join(
    sur_table %>%
      select(
        equation,
        term,
        SUR_estimate = estimate,
        SUR_SE = std_error
      ),
    by = c("equation", "term")
  ) %>%
  mutate(
    SE_change_percent = 100 * (SUR_SE / OLS_SE - 1)
  )

print(final_comparison)

# End of SUR-lecture.R
