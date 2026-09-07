# =============================================================================
# CODING SEMINAR: SEEMINGLY UNRELATED REGRESSIONS (SUR)
# =============================================================================
#
# This script consolidates and revises four earlier Chapter 11 scripts.
# It is designed for a coding seminar of approximately 90 minutes.
#
# Structure
# ---------
# Part 1: General Electric and Westinghouse investment equations
#         - separate OLS
#         - contemporaneous-correlation test
#         - system OLS versus SUR
#         - direct connection to the GLS matrix formula
#
# Part 2: A three-commodity demand system
#         - SUR with different regressors across equations
#         - the identical-regressor result
#         - cross-equation economic restrictions
#
# Part 3: Truffle demand and supply revisited
#         - equation-by-equation 2SLS
#         - system estimation by 3SLS
#         - why 3SLS is not the same estimator as SUR
#
# The script intentionally omits the panel-data, Fulton fish, Kmenta, and
# external seafood-data sections from the original files. Those sections
# either concern another topic or duplicate the same system-estimation ideas.
#
# Teaching convention
# -------------------
# Do not run the whole file without stopping. Pause at every section marked
# "PREDICT BEFORE RUNNING" and ask students to state what they expect.
# =============================================================================


# 0. PACKAGES

require(pacman) || {install.packages("pacman") ; require(pacman)}

p_load(systemfit, car, tidyverse)

###############################
# Calculate the Breusch-Pagan LM test for contemporaneous correlation.
#
# For M equations and T common observations:
#
#   lambda = T * sum_{i > j} r_ij^2
#
# Under H0 (all cross-equation covariances are zero), lambda has an
# asymptotic chi-squared distribution with M(M - 1)/2 degrees of freedom.
###############################

sur_lm_test <- function(...) {
  models <- list(...)

  if (length(models) < 2) {
    stop("Provide at least two fitted lm objects.")
  }

  observations <- vapply(models, nobs, numeric(1))

  if (length(unique(observations)) != 1) {
    stop("The seminar LM-test helper requires balanced equations.")
  }

  residual_matrix <- do.call(cbind, lapply(models, residuals))
  residual_correlation <- cor(residual_matrix)
  number_of_equations <- ncol(residual_matrix)
  number_of_observations <- nrow(residual_matrix)

  lower_triangle_correlations <- residual_correlation[
    lower.tri(residual_correlation)
  ]

  statistic <- number_of_observations *
    sum(lower_triangle_correlations^2)

  degrees_of_freedom <-
    number_of_equations * (number_of_equations - 1) / 2

  list(
    residuals = residual_matrix,
    correlation = residual_correlation,
    test = tibble(
      statistic = statistic,
      degrees_of_freedom = degrees_of_freedom,
      p_value = pchisq(
        statistic,
        df = degrees_of_freedom,
        lower.tail = FALSE
      )
    )
  )
}


# Extract coefficients and standard errors from a systemfit result.
extract_system_results <- function(model, estimator_name) {
  tibble(
    coefficient = names(coef(model)),
    estimator = estimator_name,
    estimate = as.numeric(coef(model)),
    std_error = sqrt(diag(vcov(model)))
  )
}


# Compare two system estimators coefficient by coefficient.
compare_system_results <- function(
    first_model,
    second_model,
    first_name,
    second_name) {

  first_results <- extract_system_results(first_model, first_name) |>
    select(
      coefficient,
      first_estimate = estimate,
      first_std_error = std_error
    )

  second_results <- extract_system_results(second_model, second_name) |>
    select(
      coefficient,
      second_estimate = estimate,
      second_std_error = std_error
    )

  left_join(first_results, second_results, by = "coefficient") |>
    mutate(
      estimate_change = second_estimate - first_estimate,
      se_change_percent =
        100 * (second_std_error / first_std_error - 1)
    )
}


# =============================================================================
# PART 1. GENERAL ELECTRIC AND WESTINGHOUSE INVESTMENT
# =============================================================================

# Economic model
# --------------
# Each firm's investment is modeled as a function of:
#
#   market value  = a proxy for expected future profitability
#   capital stock = the existing scale of the firm
#
# The two equations have different coefficient vectors. SUR does NOT pool the
# firms or force their coefficients to be equal.

investment <- tribble(
  ~period, ~investment_ge, ~market_value_ge, ~capital_ge,
           ~investment_w,  ~market_value_w,  ~capital_w,
         1,           33.1,           1170.6,        97.8,
                      12.93,            191.5,         1.8,
         2,           45.0,           2015.8,       104.4,
                      25.90,            516.0,         0.8,
         3,           77.2,           2803.3,       118.0,
                      35.05,            729.0,         7.4,
         4,           44.6,           2039.7,       156.2,
                      22.89,            560.4,        18.1,
         5,           48.1,           2256.2,       172.6,
                      18.84,            519.9,        23.5,
         6,           74.4,           2132.2,       186.6,
                      28.57,            628.5,        26.5,
         7,          113.0,           1834.1,       220.9,
                      48.51,            537.1,        36.2,
         8,           91.9,           1588.0,       287.8,
                      43.34,            561.2,        60.8,
         9,           61.3,           1749.4,       319.9,
                      37.02,            617.2,        84.4,
        10,           56.8,           1687.2,       321.3,
                      37.81,            626.7,        91.2,
        11,           93.6,           2007.7,       319.6,
                      39.27,            737.2,        92.4,
        12,          159.9,           2208.3,       346.0,
                      53.46,            760.5,        86.0,
        13,          147.2,           1656.7,       456.4,
                      55.56,            581.4,       111.1,
        14,          146.3,           1604.4,       543.4,
                      49.56,            662.3,       130.6,
        15,           98.3,           1431.8,       618.3,
                      32.04,            583.8,       141.8,
        16,           93.5,           1610.5,       647.4,
                      32.24,            635.2,       136.7,
        17,          135.2,           1819.4,       671.3,
                      54.38,            723.8,       129.7,
        18,          157.3,           2079.7,       726.1,
                      71.78,            864.1,       145.5,
        19,          179.5,           2371.6,       800.3,
                      90.08,           1193.5,       174.8,
        20,          189.6,           2759.9,       888.9,
                      68.60,           1188.9,       213.5
)

# Check that the data have the intended dimensions and no missing observations.
stopifnot(nrow(investment) == 20, !anyNA(investment))
glimpse(investment)


# 1.1 Separate equation-by-equation OLS ----------------------------------------

# PREDICT BEFORE RUNNING:
# 1. What signs do you expect on market value and capital?
# 2. Does separate OLS use shocks observed in the other firm's equation?

ge_ols <- lm(
  investment_ge ~ market_value_ge + capital_ge,
  data = investment
)

westinghouse_ols <- lm(
  investment_w ~ market_value_w + capital_w,
  data = investment
)

summary(ge_ols)
summary(westinghouse_ols)

# Expected coefficient estimates (approximately):
#
# GE:            -9.956 + 0.0266 * market value + 0.1517 * capital
# Westinghouse:  -0.509 + 0.0529 * market value + 0.0924 * capital
#
# Interpretation warning:
# These are conditional associations. OLS and SUR both require an exogeneity
# assumption for a causal interpretation.


# 1.2 Are the equation errors contemporaneously correlated? -------------------

# PREDICT BEFORE RUNNING:
# Both firms face the same macroeconomic and credit-market conditions. What
# sign would you expect for cor(e_GE,t, e_W,t)?

investment_lm <- sur_lm_test(ge_ols, westinghouse_ols)

investment_lm$correlation
investment_lm$test

# Expected values:
# residual correlation = 0.729
# LM statistic          = 10.628
# degrees of freedom    = 1
# p-value               = 0.0011
#
# Decision:
# Reject H0 that the cross-equation covariance is zero.
#
# What the test does NOT establish:
# - It does not prove that one firm's investment causes the other's.
# - It does not test for serial correlation within either equation.
# - It does not guarantee that feasible SUR has lower finite-sample MSE.

investment_residuals <- tibble(
  period = investment$period,
  `General Electric` = residuals(ge_ols),
  Westinghouse = residuals(westinghouse_ols)
) |>
  pivot_longer(
    cols = -period,
    names_to = "firm",
    values_to = "residual"
  )

ggplot(
  investment_residuals,
  aes(x = period, y = residual, colour = firm)
) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#2563A6", "#D97706")) +
  labs(
    title = "OLS residuals move together in many periods",
    subtitle = "Contemporaneous residual correlation is approximately 0.729",
    x = "Observation",
    y = "OLS residual",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")


# 1.3 Estimate the two equations as a system ----------------------------------

investment_equations <- list(
  GE = investment_ge ~ market_value_ge + capital_ge,
  Westinghouse = investment_w ~ market_value_w + capital_w
)

# DELIBERATE MISCONCEPTION:
# "Putting both equations in systemfit() and using OLS exploits their residual
# correlation."
#
# It does not. Because the stacked X matrix is block diagonal, system OLS is
# exactly equation-by-equation OLS.
investment_system_ols <- systemfit(
  investment_equations,
  method = "OLS",
  data = investment
)

separate_ols_coefficients <- c(
  GE = coef(ge_ols),
  Westinghouse = coef(westinghouse_ols)
)

max_system_ols_difference <- max(abs(
  coef(investment_system_ols) - separate_ols_coefficients
))

max_system_ols_difference
stopifnot(max_system_ols_difference < 1e-8)


# PREDICT BEFORE RUNNING:
# SUR uses an estimate of the cross-equation covariance matrix. Should it
# change only the standard errors, or can it also change the coefficients?
investment_sur <- systemfit(
  investment_equations,
  method = "SUR",
  data = investment
)

summary(investment_sur)

investment_comparison <- compare_system_results(
  investment_system_ols,
  investment_sur,
  first_name = "System OLS",
  second_name = "SUR"
)

investment_comparison

# Expected SUR coefficient estimates (approximately):
#
# GE:           -27.719 + 0.0383 * market value + 0.1390 * capital
# Westinghouse:  -1.252 + 0.0576 * market value + 0.0640 * capital
#
# SUR changes both coefficients and estimated precision. It is not merely a
# standard-error correction.


# 1.4 Connect systemfit() to the GLS matrix formula ----------------------------

# The stacked system has the form
#
#   y = X beta + e,
#   Var(e) = Omega = Sigma (x) I_T.
#
# We now reconstruct the SUR coefficient vector using the residual covariance
# matrix that systemfit actually used. This avoids differences caused only by
# alternative degrees-of-freedom corrections.

X_ge <- model.matrix(ge_ols)
X_w <- model.matrix(westinghouse_ols)
T_investment <- nrow(investment)

X_investment_system <- rbind(
  cbind(X_ge, matrix(0, T_investment, ncol(X_w))),
  cbind(matrix(0, T_investment, ncol(X_ge)), X_w)
)

y_investment_system <- c(
  investment$investment_ge,
  investment$investment_w
)

Sigma_used_by_systemfit <- investment_sur$residCovEst
Omega_used_by_systemfit <- kronecker(
  Sigma_used_by_systemfit,
  diag(T_investment)
)

Omega_inverse <- solve(Omega_used_by_systemfit)

beta_sur_from_gls <- solve(
  t(X_investment_system) %*%
    Omega_inverse %*%
    X_investment_system,
  t(X_investment_system) %*%
    Omega_inverse %*%
    y_investment_system
)

matrix_check <- max(abs(
  as.numeric(beta_sur_from_gls) -
    as.numeric(coef(investment_sur))
))

matrix_check
stopifnot(matrix_check < 1e-7)

# DISCUSS:
# Where does information cross from one equation to the other? It enters
# through Omega_inverse in the GLS normal equations. Merely stacking X and y
# and applying ordinary OLS cannot create those cross-equation weights.


# =============================================================================
# PART 2. THREE-COMMODITY DEMAND SYSTEM
# =============================================================================

# The data reproduce the textbook's three-commodity demand example.
# Variables:
#   p1, p2, p3 = commodity prices
#   income     = income
#   q1, q2, q3 = quantities demanded

demand_data <- tribble(
  ~p1,    ~p2,    ~p3, ~income,    ~q1,    ~q2,     ~q3,
  10.763,  4.474,  6.629, 487.648, 11.632, 13.194,  45.770,
  13.033, 10.836, 13.774, 364.877, 12.029,  2.181,  13.393,
   9.244,  5.856,  4.063, 541.037,  8.916,  5.586, 104.819,
   4.605, 14.010,  3.868, 760.343, 33.908,  5.231, 137.269,
  13.045, 11.417, 14.922, 421.746,  4.561, 10.930,  15.914,
   7.706,  8.755, 14.318, 578.214, 17.594, 11.854,  23.667,
   7.405,  7.317,  4.794, 561.734, 18.842, 17.045,  62.057,
   7.519,  6.360,  3.768, 301.470, 11.637,  2.682,  52.262,
   8.764,  4.188,  8.089, 379.636,  7.645, 13.008,  31.916,
  13.511,  1.996,  2.708, 478.855,  7.881, 19.623, 123.026,
   4.943,  7.268, 12.901, 433.741,  9.614,  6.534,  26.255,
   8.360,  5.839, 11.115, 525.702,  9.067,  9.397,  35.540,
   5.721,  5.160, 11.220, 513.067, 14.070, 13.188,  32.487,
   7.225,  9.145,  5.810, 408.666, 15.474,  3.340,  45.838,
   6.617,  5.034,  5.516, 192.061,  3.041,  4.716,  26.867,
  14.219,  5.926,  3.707, 462.621, 14.096, 17.141,  43.325,
   6.769,  8.187, 10.125, 312.659,  4.118,  4.695,  24.330,
   7.769,  7.193,  2.471, 400.848, 10.489,  7.639, 107.017,
   9.804, 13.315,  8.976, 392.215,  6.231,  9.089,  23.407,
  11.063,  6.874, 12.883, 377.724,  6.458, 10.346,  18.254,
   6.535, 15.533,  4.115, 343.552,  8.736,  3.901,  54.895,
  11.063,  4.477,  4.962, 301.599,  5.158,  4.350,  45.360,
   4.016,  9.231,  6.294, 294.112, 16.618,  7.371,  25.318,
   4.759,  5.907,  8.298, 365.032, 11.342,  6.507,  32.852,
   5.483,  7.077,  9.638, 256.125,  2.903,  3.770,  22.154,
   7.890,  9.942,  7.122, 184.798,  3.138,  1.360,  20.575,
   8.460,  7.043,  4.157, 359.084, 15.315,  6.497,  44.205,
   6.195,  4.142, 10.040, 629.378, 22.240, 10.963,  44.443,
   6.743,  3.369, 15.459, 306.527, 10.012, 10.140,  13.251,
  11.977,  4.806,  6.172, 347.488,  3.982,  8.637,  41.845
)

stopifnot(nrow(demand_data) == 30, !anyNA(demand_data))


# 2.1 Parsimonious equations with different regressors ------------------------

# Each equation contains its own price and income. Because p1, p2, and p3 are
# different variables, the three X matrices are not identical. This creates
# scope for SUR to change the coefficient estimates.

q1_ols <- lm(log(q1) ~ log(p1) + log(income), data = demand_data)
q2_ols <- lm(log(q2) ~ log(p2) + log(income), data = demand_data)
q3_ols <- lm(log(q3) ~ log(p3) + log(income), data = demand_data)


# PREDICT BEFORE RUNNING:
# Which signs should the own-price and income elasticities have?
summary(q1_ols)
summary(q2_ols)
summary(q3_ols)


# Test whether the three equations have contemporaneously correlated errors.
demand_lm <- sur_lm_test(q1_ols, q2_ols, q3_ols)

demand_lm$correlation
demand_lm$test

# Expected values:
# LM statistic       = 18.733
# degrees of freedom = 3
# p-value            = 0.00031

parsimonious_demand_equations <- list(
  eq1 = log(q1) ~ log(p1) + log(income),
  eq2 = log(q2) ~ log(p2) + log(income),
  eq3 = log(q3) ~ log(p3) + log(income)
)

parsimonious_demand_ols <- systemfit(
  parsimonious_demand_equations,
  method = "OLS",
  data = demand_data
)

parsimonious_demand_sur <- systemfit(
  parsimonious_demand_equations,
  method = "SUR",
  data = demand_data
)

summary(parsimonious_demand_sur)

parsimonious_demand_comparison <- compare_system_results(
  parsimonious_demand_ols,
  parsimonious_demand_sur,
  first_name = "System OLS",
  second_name = "SUR"
)

parsimonious_demand_comparison

# Expected SUR slopes (approximately):
# equation 1: own price = -0.909; income = 1.452
# equation 2: own price = -0.865; income = 1.137
# equation 3: own price = -0.999; income = 0.869


# 2.2 Test a cross-equation coefficient restriction ---------------------------

# Economic question:
# Are the three own-price elasticities equal?
#
# H0: beta_11 = beta_22 = beta_33
#
# This requires two independent restrictions. Symbolic restrictions are easier
# to read and less vulnerable to mistakes than manually locating columns in a
# large restriction matrix.

equal_own_price_restrictions <- c(
  "eq1_log(p1) - eq2_log(p2) = 0",
  "eq1_log(p1) - eq3_log(p3) = 0"
)

car::linearHypothesis(
  parsimonious_demand_sur,
  equal_own_price_restrictions,
  test = "Chisq"
)

restricted_own_price_sur <- systemfit(
  parsimonious_demand_equations,
  method = "SUR",
  data = demand_data,
  restrict.matrix = equal_own_price_restrictions
)

summary(restricted_own_price_sur)

# DISCUSS:
# Failure to reject H0 does not prove that the elasticities are identical. It
# means the sample does not provide strong evidence against the restrictions.
# If the restrictions are economically credible, imposing them can improve
# precision by reducing the number of free parameters.


# 2.3 Same regressors in every equation: the zero-gain result -----------------

# Now include all three prices and income in every demand equation. Each
# equation has exactly the same X matrix:
#
#   X1 = X2 = X3.
#
# PREDICT BEFORE RUNNING:
# The residuals can remain correlated. Will unrestricted SUR coefficients still
# differ from unrestricted equation-by-equation OLS coefficients?

full_demand_equations <- list(
  eq1 = log(q1) ~ log(p1) + log(p2) + log(p3) + log(income),
  eq2 = log(q2) ~ log(p1) + log(p2) + log(p3) + log(income),
  eq3 = log(q3) ~ log(p1) + log(p2) + log(p3) + log(income)
)

full_demand_ols <- systemfit(
  full_demand_equations,
  method = "OLS",
  data = demand_data
)

full_demand_sur <- systemfit(
  full_demand_equations,
  method = "SUR",
  data = demand_data
)

same_X_difference <- max(abs(
  coef(full_demand_ols) - coef(full_demand_sur)
))

same_X_difference
stopifnot(same_X_difference < 1e-7)

# Key conclusion:
# Correlated errors are not sufficient for an SUR efficiency gain. When all
# equations use the same regressors, unrestricted OLS and SUR coincide.


# 2.4 Optional extension: impose demand-theory restrictions -------------------

# This section is optional if seminar time is limited.
#
# Homogeneity of degree zero:
# In a log demand equation containing all prices and income, proportional
# changes in every price AND income should leave quantity demanded unchanged.
# The restriction therefore includes the income coefficient. The earlier
# scripts incorrectly summed only the three price coefficients.

homogeneity_restrictions <- c(
  paste0(
    "eq1_log(p1) + eq1_log(p2) + eq1_log(p3) + ",
    "eq1_log(income) = 0"
  ),
  paste0(
    "eq2_log(p1) + eq2_log(p2) + eq2_log(p3) + ",
    "eq2_log(income) = 0"
  ),
  paste0(
    "eq3_log(p1) + eq3_log(p2) + eq3_log(p3) + ",
    "eq3_log(income) = 0"
  )
)

# Parameter symmetry in this simplified log-linear specification.
# These equalities are pedagogical cross-equation restrictions. They should not
# be described as the most general Slutsky symmetry conditions without further
# assumptions about the demand system.
symmetry_restrictions <- c(
  "eq1_log(p2) - eq2_log(p1) = 0",
  "eq1_log(p3) - eq3_log(p1) = 0",
  "eq2_log(p3) - eq3_log(p2) = 0"
)

# Test the restrictions on the unrestricted SUR result.
homogeneity_test <- car::linearHypothesis(
  full_demand_sur,
  homogeneity_restrictions,
  test = "Chisq"
)

symmetry_test <- car::linearHypothesis(
  full_demand_sur,
  symmetry_restrictions,
  test = "Chisq"
)

combined_theory_restrictions <- c(
  homogeneity_restrictions,
  symmetry_restrictions
)

combined_restrictions_test <- car::linearHypothesis(
  full_demand_sur,
  combined_theory_restrictions,
  test = "Chisq"
)

homogeneity_test
symmetry_test
combined_restrictions_test

# Estimate the system subject to both sets of restrictions.
restricted_demand_sur <- systemfit(
  full_demand_equations,
  method = "SUR",
  data = demand_data,
  restrict.matrix = combined_theory_restrictions
)

summary(restricted_demand_sur)

# DISCUSS:
# Unrestricted SUR adds nothing to OLS when X1 = X2 = X3. Cross-equation
# restrictions are different: they explicitly connect coefficients across
# equations and can therefore change the estimates even in the equal-X case.


# =============================================================================
# PART 3. TRUFFLE DEMAND AND SUPPLY REVISITED: 2SLS VERSUS 3SLS
# =============================================================================

# This is the one deliberate extension beyond ordinary SUR. The truffle model
# is simultaneous: equilibrium price p is endogenous in both structural
# equations. Applying SUR directly would not solve that endogeneity problem.
#
# Demand: q = alpha_1 + alpha_2*p + alpha_3*ps + alpha_4*di + e_d
# Supply: q = beta_1  + beta_2*p  + beta_3*pf              + e_s
#
# Variables:
#   p  = price of premium truffles
#   q  = quantity of premium truffles
#   ps = price of substitute truffles
#   di = disposable income
#   pf = rental price of a truffle pig, a supply-cost variable


# 3.1 Load the textbook data robustly -----------------------------------------

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/truffles.def")

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/truffles.rdata"))

glimpse(truffles)

# 3.2 Estimate the structural system by 2SLS and 3SLS -------------------------

truffle_equations <- list(
  demand = q ~ p + ps + di,
  supply = q ~ p + pf
)

# The complete set of exogenous variables acts as the system instrument set.
truffle_instruments <- ~ ps + di + pf

# 2SLS instruments for endogenous price equation by equation. It does not use
# the cross-equation covariance matrix to estimate the coefficient vectors.
truffle_2sls <- systemfit(
  truffle_equations,
  method = "2SLS",
  inst = truffle_instruments,
  data = truffles
)

# 3SLS adds a system-GLS step to the IV procedure. In that sense, it combines
# the endogeneity correction of 2SLS with the covariance weighting of SUR.
truffle_3sls <- systemfit(
  truffle_equations,
  method = "3SLS",
  inst = truffle_instruments,
  data = truffles
)

summary(truffle_2sls)
summary(truffle_3sls)

truffle_comparison <- compare_system_results(
  truffle_2sls,
  truffle_3sls,
  first_name = "2SLS",
  second_name = "3SLS"
)

truffle_comparison

# Inspect the estimated structural-error covariance matrix used by 3SLS.
truffle_3sls$residCovEst


# PREDICT AND DISCUSS:
# 1. Why would applying ordinary SUR to these two structural equations be
#    inappropriate when equilibrium price is endogenous?
# 2. Which assumption protects 2SLS and 3SLS from simultaneity bias?
# 3. Why can 3SLS be more efficient than equation-by-equation 2SLS?
# 4. Why can a misspecified equation or invalid instrument be more damaging in
#    3SLS than in equation-by-equation 2SLS?
#
# Key distinction:
#
#   SUR  = system GLS for equations whose regressors are exogenous.
#   3SLS = system IV for equations containing endogenous regressors.
#
# Both use cross-equation covariance information, but 3SLS must also solve the
# endogeneity problem with valid instruments.


# =============================================================================
# END-OF-SEMINAR CHECK
# =============================================================================

# Students should now be able to answer:
#
# 1. Why does stacked OLS reproduce equation-by-equation OLS?
# 2. What does the off-diagonal element of Sigma represent?
# 3. Which two ingredients create an SUR efficiency gain?
# 4. Why do OLS and SUR coincide when every equation has the same X matrix?
# 5. How do cross-equation restrictions differ from covariance weighting?
# 6. Why does the truffle system require 3SLS rather than ordinary SUR?

sessionInfo()

# End of SUR-coding-seminar.R
