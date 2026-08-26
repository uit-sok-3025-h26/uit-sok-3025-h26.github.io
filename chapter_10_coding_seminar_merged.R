# ==============================================================================
# CHAPTER 10 CODING SEMINAR
# Endogenous regressors, instrumental variables, and moment-based estimation
# ==============================================================================
#
# The sequence follows the theory handout:
#   Part 1: OLS and the endogeneity problem
#   Part 2: Simple IV and the mechanics of 2SLS
#   Part 3: Multiple regression IV and diagnostic tests
#   Part 4: Exercise 10.17 -- Frisch-Waugh-Lovell and instrument strength
#   Part 5: Exercise 10.18 -- alternative instruments for education
#   Part 6: Exercise 10.20 -- measurement error in the CAPM
#
# Important:
# - The code does not clear the global environment. This avoids deleting a
#   student's existing objects without warning. Restart R before the seminar if
#   you want a completely clean session.
# - Install missing packages once with, for example:
#     install.packages(c("AER", "car", "dplyr"))
# - AER::ivreg() is used throughout so that estimation and standard errors are
#   handled consistently. A manually estimated second stage gives the same
#   coefficients, but its ordinary OLS standard errors are not correct for 2SLS.
# ==============================================================================

rm(list = ls())

# 0. Setup ---------------------------------------------------------------------

library(AER)    # ivreg()
library(car)    # linearHypothesis()
library(tidyverse)  # data manipulation

# Download a Principles of Econometrics data file and return its data frame.
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mroz.rdata"))

# A small helper for the classical, homoskedastic Sargan NR^2 calculation.
# The auxiliary model must contain every instrument: included exogenous
# regressors as well as excluded instruments.
sargan_nr2 <- function(auxiliary_model, df) {
  statistic <- nobs(auxiliary_model) * summary(auxiliary_model)$r.squared
  p_value <- pchisq(statistic, df = df, lower.tail = FALSE)

  c(statistic = statistic, df = df, p_value = p_value)
}


# 1. Mroz wage data: OLS and the endogeneity problem ---------------------------
mroz_raw <- mroz

# The wage equation is defined for women with positive observed wages.
# In this data set, lfp == 1 identifies the 428 women who worked for pay.
work <- mroz_raw |>
  filter(lfp == 1)

stopifnot(nrow(work) == 428L, all(work$wage > 0))

# Structural wage equation used in the chapter:
#
#   log(wage_i) = beta_1 + beta_2 educ_i
#                 + beta_3 exper_i + beta_4 exper_i^2 + e_i
#
# The concern is that EDUC may be positively related to unobserved ability or
# other determinants contained in e. If cov(EDUC, e) != 0, OLS is inconsistent.


# Example 10.1: multiple-regression OLS ----------------------------------------

wage_ols <- lm(
  log(wage) ~ educ + exper + I(exper^2),
  data = work
)

summary(wage_ols)
confint(wage_ols, "educ")

# The coefficient on EDUC is a semi-elasticity. The common approximation is
# 100 * beta_educ percent per additional year of education. The exact percentage
# effect is 100 * (exp(beta_educ) - 1).
educ_effect_ols_approx <- 100 * coef(wage_ols)["educ"]
educ_effect_ols_exact <- 100 * (exp(coef(wage_ols)["educ"]) - 1)

c(
  approximate_percent = educ_effect_ols_approx,
  exact_percent = educ_effect_ols_exact
)


# Example 10.2: simple OLS and simple IV ---------------------------------------

wage_simple_ols <- lm(log(wage) ~ educ, data = work)
summary(wage_simple_ols)

# Candidate instrument: mother's education.
# Relevance is partly observable: the instrument must be correlated with EDUC.
# Validity is not established by this correlation. It also requires the
# instrument to be excluded from the wage equation and uncorrelated with its
# structural error.
cor(work$educ, work$mothereduc, use = "complete.obs")

wage_simple_iv_mother <- ivreg(
  log(wage) ~ educ | mothereduc,
  data = work
)

summary(wage_simple_iv_mother, diagnostics = TRUE)
confint(wage_simple_iv_mother, "educ")

# In the just-identified simple model, the IV slope is the sample analogue of
#
#   beta_IV = cov(z, y) / cov(z, x).

iv_ratio_mother <- with(
  work,
  cov(mothereduc, log(wage)) / cov(mothereduc, educ)
)

c(
  covariance_ratio = iv_ratio_mother,
  ivreg_coefficient = coef(wage_simple_iv_mother)["educ"]
)


# Example 10.3: show the two stages manually ----------------------------------

educ_first_stage_simple <- lm(educ ~ mothereduc, data = work)
work$educ_hat_simple <- fitted(educ_first_stage_simple)

wage_manual_second_stage <- lm(
  log(wage) ~ educ_hat_simple,
  data = work
)

summary(educ_first_stage_simple)
summary(wage_manual_second_stage)

c(
  manual_2sls = coef(wage_manual_second_stage)["educ_hat_simple"],
  proper_iv = coef(wage_simple_iv_mother)["educ"]
)

# The coefficients agree. Do not use the standard error printed by
# wage_manual_second_stage. Its residual is e_i + beta_2(x_i - xhat_i), not the
# structural error e_i. Use the standard error from ivreg().
c(
  naive_manual_se = coef(summary(wage_manual_second_stage))[
    "educ_hat_simple", "Std. Error"
  ],
  proper_iv_se = coef(summary(wage_simple_iv_mother))[
    "educ", "Std. Error"
  ]
)


# Example 10.4: use two excluded instruments ----------------------------------

educ_first_stage_parents <- lm(
  educ ~ mothereduc + fathereduc,
  data = work
)

summary(educ_first_stage_parents)

# This tests relevance jointly. A large F statistic does not prove validity.
linearHypothesis(
  educ_first_stage_parents,
  c("mothereduc = 0", "fathereduc = 0")
)

wage_simple_iv_parents <- ivreg(
  log(wage) ~ educ | mothereduc + fathereduc,
  data = work
)

summary(wage_simple_iv_parents, diagnostics = TRUE)


# Example 10.5: multiple-regression IV/2SLS -----------------------------------

# EDUC is treated as endogenous. EXPER and EXPER^2 are treated as exogenous and
# must appear on both sides of the vertical bar. Included exogenous variables
# act as their own instruments.
wage_iv_parents <- ivreg(
  log(wage) ~ educ + exper + I(exper^2) |
    exper + I(exper^2) + mothereduc + fathereduc,
  data = work
)

summary(wage_iv_parents, diagnostics = TRUE)
confint(wage_iv_parents, "educ")

# Compare OLS and IV estimates, but remember the trade-off:
# - OLS is more precise if EDUC is exogenous;
# - OLS is inconsistent if EDUC is endogenous;
# - IV is consistent only if the instruments are valid and relevant;
# - IV commonly has a larger standard error.
rbind(
  OLS = c(
    estimate = coef(wage_ols)["educ"],
    std_error = coef(summary(wage_ols))["educ", "Std. Error"]
  ),
  IV_2SLS = c(
    estimate = coef(wage_iv_parents)["educ"],
    std_error = coef(summary(wage_iv_parents))["educ", "Std. Error"]
  )
)


# Example 10.7: regression-based Hausman test ---------------------------------

# First stage for the potentially endogenous variable EDUC.
educ_first_stage <- lm(
  educ ~ exper + I(exper^2) + mothereduc + fathereduc,
  data = work
)

work$vhat_educ <- residuals(educ_first_stage)

# Artificial/control-function regression:
#
#   log(wage) = beta_1 + beta_2 educ + controls + delta vhat + error.
#
# H0: delta = 0  -> EDUC is exogenous.
# H1: delta != 0 -> EDUC is endogenous, provided the instruments are valid.
wage_hausman <- lm(
  log(wage) ~ educ + exper + I(exper^2) + vhat_educ,
  data = work
)

summary(wage_hausman)
linearHypothesis(wage_hausman, "vhat_educ = 0")

# A failure to reject does not prove exogeneity. In particular, the test can
# have low power when the instruments are weak.


# Example 10.7: test the surplus moment restriction ----------------------------

# There are L = 2 excluded instruments and B = 1 endogenous regressor, so the
# overidentification test has L - B = 1 degree of freedom.
work$ehat_iv_parents <- residuals(wage_iv_parents)

sargan_aux_mroz <- lm(
  ehat_iv_parents ~ exper + I(exper^2) + mothereduc + fathereduc,
  data = work
)

sargan_mroz <- sargan_nr2(sargan_aux_mroz, df = 2 - 1)
sargan_mroz

# Under homoskedasticity, N*R^2 is asymptotically chi-squared. A large statistic
# rejects the joint validity of the surplus restrictions. It does not tell us
# which instrument is invalid. The diagnostics in ivreg() also report this test.
summary(wage_iv_parents, diagnostics = TRUE)


# 2. Exercise 10.17: FWL and first-stage strength ------------------------------

# This exercise makes equation (10.25) concrete. It uses only MOTHEREDUC as the
# excluded instrument and controls for EXPER and EXPER^2.

# (a) Full first stage. REDUCHAT is the residual from the unrestricted model.
educ_first_stage_mother <- lm(
  educ ~ exper + I(exper^2) + mothereduc,
  data = work
)

work$REDUCHAT <- residuals(educ_first_stage_mother)
sse_full_first_stage <- deviance(educ_first_stage_mother)
sse_full_first_stage

# (b) Residualize EDUC with respect to the included exogenous controls.
educ_on_controls <- lm(educ ~ exper + I(exper^2), data = work)
work$REDUC <- residuals(educ_on_controls)
sse_educ_controls <- deviance(educ_on_controls)
sse_educ_controls

# (c) Residualize MOTHEREDUC with respect to the same controls.
mothereduc_on_controls <- lm(
  mothereduc ~ exper + I(exper^2),
  data = work
)

work$RMOM <- residuals(mothereduc_on_controls)
ss_rmom <- sum(work$RMOM^2)
ss_rmom

# (d) Frisch-Waugh-Lovell regression. Because both variables are residuals from
# regressions containing an intercept, they have means of zero. The no-intercept
# form makes the intended geometry explicit.
fwl_first_stage <- lm(REDUC ~ 0 + RMOM, data = work)
summary(fwl_first_stage)

# FWL identities: the partial-regression slope equals the coefficient on
# MOTHEREDUC in the full first stage, and the residual sums of squares agree.
c(
  full_first_stage = coef(educ_first_stage_mother)["mothereduc"],
  fwl_regression = coef(fwl_first_stage)["RMOM"]
)

c(
  full_first_stage_sse = sse_full_first_stage,
  fwl_sse = deviance(fwl_first_stage)
)

# The variation in EDUC explained uniquely by the excluded instrument is
#
#   theta_hat^2 * sum(RMOM_i^2).
#
# This is the denominator component in the simple form of equation (10.25).
theta_hat <- coef(fwl_first_stage)["RMOM"]
instrument_induced_variation <- theta_hat^2 * ss_rmom
partial_r_squared <- instrument_induced_variation / sse_educ_controls

c(
  theta_squared_times_ss_rmom = instrument_induced_variation,
  reduction_in_sse = sse_educ_controls - sse_full_first_stage,
  partial_r_squared = partial_r_squared
)

# Under the chapter's homoskedastic assumptions, the approximate variance of
# the IV coefficient on EDUC is inversely related to this instrument-induced
# variation:
#
#   Var(beta_hat_educ^IV) = sigma^2 /
#                           [theta_hat^2 * sum(RMOM_i^2)].
#
# A small theta_hat or little residual variation in MOTHEREDUC makes the IV
# estimate imprecise. This is the mathematical meaning of a weak instrument.


# 3. Exercise 10.18: parental college indicators as instruments ----------------

work <- work |>
  mutate(
    mothercoll = as.integer(mothereduc > 12),
    fathercoll = as.integer(fathereduc > 12)
  )

# (a) Percentage of parents with some college education.
work |>
  summarise(
    pct_mothercoll = 100 * mean(mothercoll),
    pct_fathercoll = 100 * mean(fathercoll)
  )

# (b) Pairwise correlations. Correlation with EDUC concerns relevance only.
cor(
  work[, c("educ", "mothercoll", "fathercoll")],
  use = "complete.obs"
)

# A larger raw correlation does not automatically make an instrument better.
# Strength is conditional on the controls, and instrument choice must also be
# supported by credible exogeneity and exclusion arguments.

# (c) One excluded instrument: MOTHERCOLL. The model is just identified.
wage_iv_mothercoll <- ivreg(
  log(wage) ~ educ + exper + I(exper^2) |
    exper + I(exper^2) + mothercoll,
  data = work
)

summary(wage_iv_mothercoll, diagnostics = TRUE)
confint(wage_iv_mothercoll, "educ")

# (d) First-stage strength of MOTHERCOLL conditional on the controls.
educ_first_stage_mothercoll <- lm(
  educ ~ exper + I(exper^2) + mothercoll,
  data = work
)

summary(educ_first_stage_mothercoll)
linearHypothesis(educ_first_stage_mothercoll, "mothercoll = 0")

# With one excluded instrument, t^2 for its coefficient equals the one-variable
# partial F statistic. The F > 10 rule is only a rough diagnostic, not a proof
# that weak identification is absent.

# (e) Two excluded instruments: MOTHERCOLL and FATHERCOLL.
wage_iv_parentcoll <- ivreg(
  log(wage) ~ educ + exper + I(exper^2) |
    exper + I(exper^2) + mothercoll + fathercoll,
  data = work
)

summary(wage_iv_parentcoll, diagnostics = TRUE)
confint(wage_iv_parentcoll, "educ")

# Compare the confidence intervals using one and two excluded instruments.
rbind(
  MOTHERCOLL = confint(wage_iv_mothercoll, "educ"),
  BOTH_PARENT_INDICATORS = confint(wage_iv_parentcoll, "educ")
)

# (f) Joint first-stage significance of both excluded instruments.
educ_first_stage_parentcoll <- lm(
  educ ~ exper + I(exper^2) + mothercoll + fathercoll,
  data = work
)

summary(educ_first_stage_parentcoll)
linearHypothesis(
  educ_first_stage_parentcoll,
  c("mothercoll = 0", "fathercoll = 0")
)

# (g) Manual Sargan test. Here L - B = 2 - 1 = 1.
work$ehat_iv_parentcoll <- residuals(wage_iv_parentcoll)

sargan_aux_parentcoll <- lm(
  ehat_iv_parentcoll ~ exper + I(exper^2) + mothercoll + fathercoll,
  data = work
)

sargan_parentcoll <- sargan_nr2(sargan_aux_parentcoll, df = 1)
sargan_parentcoll

# This is a joint test of the surplus moment restriction. Failure to reject is
# not proof that both instruments are valid; the test can have limited power.


# 4. Exercise 10.20: CAPM and measurement error -------------------------------
#browseURL("https://www.principlesofeconometrics.com/poe5/data/def/capm5.def")

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/capm5.rdata"))

capm <- capm5 |>
  mutate(
    market_premium = mkt - riskfree,
    msft_premium = msft - riskfree
  )

# CAPM equation:
#
#   r_msft - r_f = alpha + beta(r_m - r_f) + e.
#
# If the observed market premium measures the true market portfolio with
# classical error, OLS beta is attenuated toward zero in the simple-regression
# model. The proposed constructed instruments below are useful pedagogically,
# but their validity should be treated with great caution.


# (a) Estimate Microsoft beta by OLS -------------------------------------------

capm_ols <- lm(msft_premium ~ market_premium, data = capm)
summary(capm_ols)
confint(capm_ols, "market_premium")

# beta > 1 indicates that the security tends to move more than the measured
# market portfolio; beta < 1 indicates a more defensive security.
linearHypothesis(capm_ols, "market_premium = 1")


# (b) Construct RANK and assess first-stage relevance --------------------------

# rank() avoids changing the chronological row order of the data.
capm <- capm |>
  mutate(rank_market = rank(market_premium, ties.method = "average"))

market_first_stage_rank <- lm(
  market_premium ~ rank_market,
  data = capm
)

summary(market_first_stage_rank)
linearHypothesis(market_first_stage_rank, "rank_market = 0")

# RANK is strongly related to the market premium by construction. This almost
# guarantees first-stage relevance, but it does NOT establish validity. Because
# RANK is a deterministic transformation of the potentially mismeasured
# endogenous regressor, it may inherit its correlation with the structural
# error. A huge first-stage F statistic cannot repair this problem.


# (c) Regression-based Hausman test using RANK --------------------------------

capm$vhat_rank <- residuals(market_first_stage_rank)

capm_hausman_rank <- lm(
  msft_premium ~ market_premium + vhat_rank,
  data = capm
)

summary(capm_hausman_rank)
linearHypothesis(capm_hausman_rank, "vhat_rank = 0")

# The exercise asks for a 1% test. Compare the p-value with 0.01. The conclusion
# is conditional on RANK being a valid instrument—a demanding assumption here.


# (d) IV/2SLS using RANK -------------------------------------------------------

capm_iv_rank <- ivreg(
  msft_premium ~ market_premium | rank_market,
  data = capm
)

summary(capm_iv_rank, diagnostics = TRUE)
confint(capm_iv_rank, "market_premium")

rbind(
  OLS = c(
    estimate = coef(capm_ols)["market_premium"],
    std_error = coef(summary(capm_ols))["market_premium", "Std. Error"]
  ),
  IV_RANK = c(
    estimate = coef(capm_iv_rank)["market_premium"],
    std_error = coef(summary(capm_iv_rank))[
      "market_premium", "Std. Error"
    ]
  )
)

# Under classical measurement error, one expects OLS attenuation and potentially
# a larger IV coefficient. That qualitative agreement is not evidence that the
# instrument is valid.


# (e) Add POS as a second constructed instrument -------------------------------

capm <- capm |>
  mutate(pos_market = as.integer(market_premium > 0))

market_first_stage_both <- lm(
  market_premium ~ rank_market + pos_market,
  data = capm
)

summary(market_first_stage_both)
linearHypothesis(
  market_first_stage_both,
  c("rank_market = 0", "pos_market = 0")
)

# The joint F statistic assesses relevance. Like RANK, POS is constructed from
# the endogenous regressor itself, so its exogeneity is not guaranteed.


# (f) Hausman test using the two-instrument first stage -------------------------

capm$vhat_both <- residuals(market_first_stage_both)

capm_hausman_both <- lm(
  msft_premium ~ market_premium + vhat_both,
  data = capm
)

summary(capm_hausman_both)
linearHypothesis(capm_hausman_both, "vhat_both = 0")


# (g) IV/2SLS using RANK and POS -----------------------------------------------

capm_iv_both <- ivreg(
  msft_premium ~ market_premium | rank_market + pos_market,
  data = capm
)

summary(capm_iv_both, diagnostics = TRUE)
confint(capm_iv_both, "market_premium")

rbind(
  OLS = c(
    estimate = coef(capm_ols)["market_premium"],
    std_error = coef(summary(capm_ols))["market_premium", "Std. Error"]
  ),
  IV_RANK_POS = c(
    estimate = coef(capm_iv_both)["market_premium"],
    std_error = coef(summary(capm_iv_both))[
      "market_premium", "Std. Error"
    ]
  )
)


# (h) Manual Sargan test of the surplus instrument -----------------------------

# L = 2 instruments and B = 1 endogenous regressor, so df = L - B = 1.
capm$ehat_iv_both <- residuals(capm_iv_both)

sargan_aux_capm <- lm(
  ehat_iv_both ~ rank_market + pos_market,
  data = capm
)

sargan_capm <- sargan_nr2(sargan_aux_capm, df = 1)
sargan_capm
qchisq(0.95, df = 1)

# A rejection means that at least one surplus moment restriction is inconsistent
# with the data under the maintained model. It does not identify the offending
# instrument. Here, the rejection reinforces the conceptual concern that RANK
# and POS are functions of the potentially endogenous measured market return.


# 5. Closing checklist ----------------------------------------------------------

# Before interpreting an IV coefficient causally, answer all four questions:
#
# 1. Endogeneity:
#    Why is an included regressor correlated with the structural error?
#
# 2. Relevance:
#    Do the excluded instruments explain the endogenous regressor conditional on
#    all included exogenous controls? Inspect coefficients, the partial F, and
#    the partial R-squared—not only raw correlations.
#
# 3. Exogeneity and exclusion:
#    Why should each instrument be uncorrelated with the structural error and
#    affect the outcome only through the endogenous regressor? These are mainly
#    economic or institutional arguments, not conclusions from a first stage.
#
# 4. Diagnostics and uncertainty:
#    Are the instruments sufficiently strong? Is the model overidentified, and
#    if so, what does the overidentification test show? Are robust, clustered, or
#    serial-correlation-robust standard errors needed for the application?
#
# The central lesson of Chapter 10 is that IV replaces the OLS moment condition
# E[x_i e_i] = 0 with instrument moment conditions E[z_i e_i] = 0. The code can
# estimate those moments, but the credibility of the instrument must come from
# the research design.
