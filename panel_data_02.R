# =============================================================================
# panel_data_02.R
# Random Effects and Model Choice · Which variation can we trust?
#
# Accompanies panel_data_02.qmd and follows the style of panel_data_01.R.
# Hill, Griffiths & Lim, Principles of Econometrics (5th ed.), Chapter 15,
# §15.3 and §§15.4–15.4.3, 15.4.5; handout Version 2.1.
#
# Lab sequence: 2.1 clustering · 2.2 LM test · 2.3 RE weighting ·
# 2.4 conventional Hausman test · 2.5 Mundlak decomposition ·
# 2.6 cluster-robust Mundlak test · 2.7 reporting.
#
# Run one numbered section at a time. Before each lab, make a prediction.
# Blocks marked [EXTRA] are optional reference material for self-study.
# This script runs independently of Lecture 1. Packages must be installed.
# =============================================================================

rm(list = ls())  # clear workspace
getwd()

# Packages --------------------------------------------------------------------
library(tidyverse)
library(plm)       # panel estimators; masks dplyr::between, lag and lead
library(sandwich)  # cluster-robust covariance matrices
library(lmtest)    # coeftest() and coefci()

# Data ------------------------------------------------------------------------
# Definition: https://www.principlesofeconometrics.com/poe5/data/def/liquor5.def
# hh = household; year = period coded 1–3 (not calendar years).
# liquor and income are both measured in thousands of dollars.
load(url("https://www.principlesofeconometrics.com/poe5/data/rdata/liquor5.rdata"))

panel <- as_tibble(liquor5) |>
  arrange(hh, year) |>
  group_by(hh) |>
  mutate(
    income_mean   = mean(income),
    liquor_mean   = mean(liquor),
    income_within = income - income_mean,
    liquor_within = liquor - liquor_mean
  ) |>
  ungroup()

panel
count(panel, hh) |> count(n, name = "households")
N  <- n_distinct(panel$hh)
TT <- n_distinct(panel$year)  # TT, because T is shorthand for TRUE
NT <- nrow(panel)
G  <- N                     # households are the clusters
c(N = N, T = TT, NT = NT)     # 40 households, 3 periods, 120 observations

# 0. Retrieval and connection -------------------------------------------------
# Why might the within estimate be imprecise with 120 observations?
Wxx <- sum(panel$income_within^2)
Wxy <- sum(panel$income_within * panel$liquor_within)
Txx <- sum((panel$income - mean(panel$income))^2)
Txy <- sum((panel$income - mean(panel$income)) *
           (panel$liquor - mean(panel$liquor)))
Bxx <- Txx - Wxx
Bxy <- Txy - Wxy
within_share <- Wxx / Txx
within_share  # about 0.0378: only 3.78% of income variation is within households

household_means <- panel |> distinct(hh, income_mean, liquor_mean)
between <- lm(liquor_mean ~ income_mean, data = household_means)
within_manual <- lm(liquor_within ~ 0 + income_within, data = panel)
b_within  <- coef(within_manual)[["income_within"]]
b_between <- coef(between)[["income_mean"]]

# 1. The hook -----------------------------------------------------------------
# The lecture reports FE = 0.02074 (clustered SE 0.01943) and
# RE = 0.02658 (clustered SE 0.00746). We reproduce these below.
# Is a smaller SE sufficient reason to use RE? State its extra assumption:
# E(u_i | income_i1, ..., income_iT) = 0.
# Strict exogeneity of the time-varying error is still required.

# 3. The econometric lab ------------------------------------------------------

# Lab 2.1 · What does clustering change? --------------------------------------
# Prediction: does allowing dependence within households change the slope?
pooled <- lm(liquor ~ income, data = panel)
summary(pooled)$coefficients

# Same convention as Lecture 1: CR1 = G/(G-1) * (NT-1)/(NT-k),
# where k counts the columns of the fitted regression's design matrix.
# vcovCL(type = "HC1", cadjust = TRUE) applies both factors.
# Use t(G-1) = t(39) for tests and confidence intervals throughout.
vc_pooled_cluster <- vcovCL(pooled, cluster = panel$hh,
                            type = "HC1", cadjust = TRUE)
pooled_inference <- coeftest(pooled, vcov. = vc_pooled_cluster, df = G - 1)
pooled_inference
coefci(pooled, vcov. = vc_pooled_cluster, df = G - 1)

# Slope stays 0.027073; SE changes from 0.00533 to 0.007684.
# Clustering changes uncertainty, not the estimated coefficient.
# It does NOT correct correlation between income and omitted preferences.
# Clustered SEs need not be larger in every sample.

# Lab 2.2 · Is there a household component? ----------------------------------
# Prediction: clustering raised the pooled SE in Lab 2.1. What does that
# suggest about the residual cross-products in POE5 eq. (15.35)?
# H0: sigma_u2 = 0 (pooled model adequate)  vs  H1: sigma_u2 > 0.
# The LM test needs only the pooled OLS residuals.
panel <- panel |>
  mutate(e_pooled = resid(pooled))

household_sums <- panel |>
  group_by(hh) |>
  summarise(e_sum = sum(e_pooled), .groups = "drop")

# Squaring each household's residual sum adds cross-products e_is * e_it.
# They are near zero without a household component and positive with one.
ratio <- sum(household_sums$e_sum^2) / sum(panel$e_pooled^2)
LM <- sqrt(NT / (2 * (TT - 1))) * (ratio - 1)
c(ratio = ratio, LM = LM, p_one_sided = pnorm(LM, lower.tail = FALSE))
# LM is about 4.55 > 1.645: reject sigma_u2 = 0 at the 5% level.
# Under H0, LM ~ N(0, 1); the test is ONE-sided. A negative LM is not
# evidence of household effects.

# Cross-check with plm.
pooled_panel <- plm(liquor ~ income, data = panel, index = c("hh", "year"),
                    model = "pooling")
plmtest(pooled_panel, effect = "individual", type = "honda")  # LM, one-sided
plmtest(pooled_panel, effect = "individual", type = "bp")     # LM^2, chi2(1)
LM^2  # the "bp" statistic
# POE5 footnote 11: with LM^2, the one-sided 5% critical value is
# 2.706 = 1.645^2, and it applies only when LM > 0.

# Interpretation: residuals from the same household move together. This is
# the dependence that made the conventional pooled SE too small in Lab 2.1.
# The LM test does NOT say whether household effects are related to income,
# so it cannot choose between RE and FE. That is the job of Labs 2.4 and 2.6.

# Lab 2.3 · How far does RE move towards FE? ----------------------------------
# Prediction: does downweighting between variation put RE close to FE?
random <- plm(liquor ~ income, data = panel, index = c("hh", "year"),
              model = "random", effect = "individual", random.method = "swar")
summary(random)
# summary(random) reports conventional RE inference. Clustered inference
# for the lecture's comparison is calculated below on quasi-demeaned data.

# Obtain the estimated variance components directly, as in the handout.
sigma_e2 <- random$ercomp$sigma2[["idios"]]
sigma_u2 <- random$ercomp$sigma2[["id"]]
omega  <- sigma_e2 / (sigma_e2 + TT * sigma_u2)
lambda <- 1 - sqrt(omega)
c(sigma_e2 = sigma_e2, sigma_u2 = sigma_u2, omega = omega, lambda = lambda)
# Approximately 0.96403, 0.72512, 0.30708 and 0.44586.
# Notation bridge to handout V2.1:
# sigma_u2 = sigma_alpha^2; omega = phi^2; lambda = theta (plm's theta).
random$ercomp$theta

# Quasi-demean: remove a fraction lambda of each household mean.
# Transform the intercept too: its regressor becomes 1 - lambda.
panel <- panel |>
  mutate(
    constant_star = 1 - lambda,
    income_star   = income - lambda * income_mean,
    liquor_star   = liquor - lambda * liquor_mean
  )
random_manual <- lm(liquor_star ~ 0 + constant_star + income_star, data = panel)
coef(random_manual)
coef(random)
# The constant_star coefficient is the intercept in the original RE model.
# Dropping constant_star and fitting through zero generally changes the slope.

# Same slope from within and between sums of squares (balanced, one regressor).
b_pooled <- coef(pooled)[["income"]]
b_random <- coef(random_manual)[["income_star"]]
b_re_weighted <- (Wxy + omega * Bxy) / (Wxx + omega * Bxx)
c(quasi_demeaning = b_random, sums_of_squares = b_re_weighted,
  plm = coef(random)[["income"]])  # all about 0.026575

# omega is NOT the final coefficient weight on the between slope.
weight_between <- omega * Bxx / (Wxx + omega * Bxx)
weight_within  <- 1 - weight_between
c(within_weight = weight_within, between_weight = weight_between)
# Between still receives about 88.7% of the coefficient weight!
weight_within * b_within + weight_between * b_between

# Figure: same sensitivity relationship as the lecture's slider.
# We vary omega holding the sample sums of squares fixed. These are not
# separately estimated models, and the curve is not an exogeneity test.
weight_curve <- tibble(weight = seq(0, 1, by = 0.001)) |>
  mutate(slope = (Wxy + weight * Bxy) / (Wxx + weight * Bxx))

fig_re_weight <- ggplot(weight_curve, aes(x = weight, y = slope)) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  geom_hline(yintercept = b_between, linetype = "dashed", colour = "grey45") +
  annotate("point", x = c(0, omega, 1),
           y = c(b_within, b_random, b_pooled), size = 2.5) +
  annotate("text", x = c(0, omega, 1),
           y = c(b_within, b_random, b_pooled),
           label = c("FE", "RE", "Pooled OLS"),
           hjust = c(-0.3, -0.1, 1.1), vjust = c(0.5, 1.8, 1.8)) +
  annotate("text", x = 0.02, y = b_between, label = "Between",
           hjust = 0, vjust = -0.6) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = expression(omega ~ "(between sums-of-squares weight)"),
       y = "Income slope", title = "How much does RE move towards FE?") +
  theme_minimal(base_size = 12)
fig_re_weight

# Interpretation: a large Bxx keeps RE close to pooling despite downweighting.
# Similar pooled and RE estimates do not independently confirm causality.
# This scalar weighting identity does not imply coefficient ordering in
# general multiple-regressor models.

# [EXTRA] Where do the Swamy–Arora variance components come from?
# In this balanced example, the within residual variance uses NT - N - 1 df.
# The between residual variance estimates sigma_u2 + sigma_e2/TT.
sigma_e2_manual <- sum(resid(within_manual)^2) / (NT - N - 1)
sigma_u2_manual <- max(0, sum(resid(between)^2) / (N - 2) - sigma_e2_manual / TT)
c(within_residual_variance = sigma_e2_manual,
  household_variance = sigma_u2_manual)
# max(0, ...) imposes a nonnegative household variance; it is not binding here.

# Lab 2.4 · Do FE and RE differ? The conventional Hausman test ----------------
# Prediction: the slopes differ by about 0.006. Is that large?
# H0: household effects unrelated to income -> FE and RE both consistent,
#     and RE is efficient.
# H1: only FE is consistent.
fixed <- plm(liquor ~ income, data = panel, index = c("hh", "year"),
             model = "within", effect = "individual")

b_fe  <- coef(fixed)[["income"]]
b_re  <- coef(random)[["income"]]
se_fe <- sqrt(vcov(fixed)[["income", "income"]])    # conventional SE
se_re <- sqrt(vcov(random)[["income", "income"]])   # conventional SE

# POE5 eq. (15.36). The denominator uses eq. (15.37):
# var(b_FE - b_RE) = var(b_FE) - var(b_RE),
# which holds because cov(b_FE, b_RE) = var(b_RE) when RE is EFFICIENT.
var_difference <- se_fe^2 - se_re^2
var_difference  # can be negative in finite samples; then the test fails
t_hausman <- (b_fe - b_re) / sqrt(var_difference)
c(b_fe = b_fe, b_re = b_re, se_fe = se_fe, se_re = se_re,
  t = t_hausman, p = 2 * pnorm(-abs(t_hausman)))
# About: SEs 0.0209 and 0.0070, t about -0.30, p about 0.77.

phtest(fixed, random)
t_hausman^2  # with one regressor, the chi2(1) statistic equals t^2

# Do not reject under the conventional assumptions (homoskedastic,
# serially uncorrelated e_it).

# Deliberate misconception: insert the CLUSTERED SEs into eq. (15.36).
t_wrong <- (0.02074 - 0.02658) / sqrt(0.01943^2 - 0.00746^2)
t_wrong  # about -0.33: looks similar, but has no justification.
# With clustered SEs we no longer assume RE is efficient, so eq. (15.37)
# fails and the denominator is not the variance of the difference.
# The cluster-robust alternative is the Mundlak test in Lab 2.6.

# Lab 2.5 · What happens when mean income is included? -------------------------
# Prediction: which coefficient should reproduce FE?
mundlak <- lm(liquor ~ income + income_mean, data = panel)
vc_mundlak_cluster <- vcovCL(mundlak, cluster = panel$hh,
                             type = "HC1", cadjust = TRUE)
mundlak_inference <- coeftest(mundlak, vcov. = vc_mundlak_cluster, df = G - 1)
mundlak_inference

# y_it = delta + beta * income_it + gamma * mean_income_i + error_it
#      = delta + beta * income_within_it
#              + (beta + gamma) * mean_income_i + error_it.
c(within = coef(mundlak)[["income"]],
  between = coef(mundlak)[["income"]] + coef(mundlak)[["income_mean"]],
  FE = b_within, between_regression = b_between)
# beta = 0.020742; gamma = 0.006579; beta + gamma = 0.027321.
# This augmented OLS specification illustrates the Mundlak decomposition.
# gamma = between slope minus within slope; gamma = 0 means the two sources
# of variation tell the same story.
# The income coefficient equals FE; its CR1 SE differs slightly because
# k = 3 here, compared with k = 1 in the demeaned FE slope regression.

# [EXTRA] RE on the augmented equation gives the same beta and gamma (POE5
# Tables 15.6-15.7). Quasi-demean the mean regressor as well: it becomes
# (1 - lambda) * income_mean. This avoids plm's between regression, where
# income and income_mean would be perfectly collinear.
panel <- panel |>
  mutate(income_mean_star = (1 - lambda) * income_mean)
mundlak_re <- lm(liquor_star ~ 0 + constant_star + income_star +
                   income_mean_star, data = panel)
cbind(OLS = coef(mundlak),
      RE  = coef(mundlak_re)[c("constant_star", "income_star",
                               "income_mean_star")])
# Identical for any lambda: the intercept and the mean regressor are scaled
# by the same factor, so beta is still identified only from within changes.
# Mundlak is also a MODEL (POE5 §15.4.6): with the means included, RE can
# keep time-invariant regressors that FE drops. There are none in liquor5.

# [EXTRA] Explicit within–between parameterisation: same fitted values.
within_between <- lm(liquor ~ income_within + income_mean, data = panel)
coef(within_between)
# Here the income_mean coefficient is the BETWEEN slope, not gamma.

# Lab 2.6 · Does the robust diagnostic settle model choice? -------------------
# H0: gamma = 0. Use the clustered t-test on the mean-income coefficient.
mundlak_test <- mundlak_inference["income_mean", , drop = FALSE]
mundlak_test
coefci(mundlak, vcov. = vc_mundlak_cluster, df = G - 1)["income_mean", ]
# Estimate 0.006579; clustered SE 0.020807; t(39); p about 0.7535.
# We do not reject this restriction in the specified Mundlak representation.
# It agrees with the conventional Hausman test (Lab 2.4, p about 0.77), but
# only this version is valid without homoskedastic, uncorrelated e_it.
c(hausman_conventional = 2 * pnorm(-abs(t_hausman)),
  mundlak_clustered = mundlak_test[, "Pr(>|t|)"])
# With several regressors: add one mean per time-varying regressor and use
# a joint Wald test, e.g. lmtest::waldtest(..., vcov = ...), chi2(K_S).
# Non-rejection neither proves independence nor establishes strict exogeneity.
# The within–between difference is uncertain; substantive reasoning matters.

# Deliberate misconception: "p = 0.75 proves RE is valid."
# Replace it with a statement that distinguishes lack of evidence from proof.

# Lab 2.7 · What would you report? --------------------------------------------
# Prediction: choose the economic interpretation before comparing intervals.
vc_fe_cluster <- vcovCL(within_manual, cluster = panel$hh,
                        type = "HC1", cadjust = TRUE)
vc_re_cluster <- vcovCL(random_manual, cluster = panel$hh,
                        type = "HC1", cadjust = TRUE)
fe_inference <- coeftest(within_manual, vcov. = vc_fe_cluster, df = G - 1)
re_inference <- coeftest(random_manual, vcov. = vc_re_cluster, df = G - 1)
fe_inference
re_inference

# k = 2 for pooled and quasi-demeaned RE; k = 1 for demeaned FE.
# Cluster the transformed RE regression to reproduce the lecture convention.
# As usual for feasible GLS, the fitted transformation is treated as given.
comparison <- tibble(
  model = c("Pooled", "FE", "RE"),
  estimate = c(b_pooled, b_within, b_random),
  std_error = c(sqrt(vc_pooled_cluster["income", "income"]),
                sqrt(vc_fe_cluster["income_within", "income_within"]),
                sqrt(vc_re_cluster["income_star", "income_star"]))
) |>
  mutate(
    statistic = estimate / std_error,
    p_value = 2 * pt(-abs(statistic), df = G - 1),
    conf_low = estimate - qt(0.975, df = G - 1) * std_error,
    conf_high = estimate + qt(0.975, df = G - 1) * std_error
  )
comparison |> mutate(across(where(is.numeric), ~ round(.x, 6)))
#             estimate       SE         95% interval
# Pooled      0.027073     0.007684     [ 0.011530, 0.042615]
# FE          0.020742     0.019427     [-0.018553, 0.060037]
# RE          0.026575     0.007458     [ 0.011489, 0.041662]

fig_coefficient_intervals <- comparison |>
  mutate(model = factor(model, levels = c("RE", "FE", "Pooled"))) |>
  ggplot(aes(x = model, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high), colour = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Income slope and 95% confidence interval",
       title = "Which variation can we trust?",
       caption = "Household-clustered CR1; t(39) intervals.") +
  theme_minimal(base_size = 12)
fig_coefficient_intervals

# Both variables are in $1,000. Multiply slopes and interval endpoints by
# 1,000 to express dollars of liquor expenditure per $1,000 more income.
comparison |>
  select(model, estimate, conf_low, conf_high) |>
  mutate(across(where(is.numeric), ~ 1000 * .x))

# Report a positive but uncertain within-household association. RE gains
# precision by using between information under stronger restrictions.
# Neither a smaller interval nor a large diagnostic p-value proves causality.
# Forty clusters support approximate inference, not an exact finite-sample test.

# Save figures ----------------------------------------------------------------

# Install if you don't have it: install.packages("this.path")
library(this.path)

# Sets directory to the folder containing the script
setwd(this.dir())
getwd()

figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "re_weight.png"), fig_re_weight,
       width = 7, height = 4, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "coefficient_intervals.png"), fig_coefficient_intervals,
       width = 7, height = 4, dpi = 300, bg = "white")

# 4. The AI audit -------------------------------------------------------------
# Use these results in the lecture's audit prompt:
LM
comparison
mundlak_test
# Evaluate: "The LM test shows random effects are present, and the large
# Mundlak p-value proves RE is valid. Because RE is more precise, we should
# report only RE and interpret its coefficient causally."
# Check the answer separates the LM result (presence of u_i), the RE
# restriction (relation to income), precision and causality.
# Stress test: should a weak diagnostic overrule credible evidence that
# stable preferences and income are related? Explain why not.

# 5. Exit problem -------------------------------------------------------------
# "Clustering RE errors by firm makes RE consistent even when managerial
# ability is correlated with employment." Correct?
# No: clustering changes uncertainty, not correlation with omitted ability.
# FE can remove stable additive ability, but time-varying confounding remains.

# [EXTRA] Mundlak with panelr · connecting the algebra to a package ------------
# Run this extension after Labs 2.5–2.7. It requires the panelr package.
# Documentation: https://panelr.jacob-long.com/reference/wbm
library(panelr)

# panel_data() records the household and time identifiers. Use the original
# variables in the formula: wbm() constructs the means/deviations for us.
liquor_panel <- panelr::panel_data(panel, id = hh, wave = year)

# A. The contextual specification: the same regressors as our Mundlak lm().
#
# liquor_it = delta + beta * income_it + gamma * mean_income_i + c_i + e_it
#
# Mundlak models the part of the household effect associated with income
# through its household mean. c_i is the remaining household random effect;
# the specification assumes its conditional mean is zero given the regressors.
# This allows the original household effect to be related to income through
# mean_income_i. It does not address time-varying endogeneity automatically.
mundlak_panelr <- panelr::wbm(
  liquor ~ income,
  data = liquor_panel,
  model = "contextual"
)
summary(mundlak_panelr)

# Read the output carefully:
# income        -> beta: the within-household slope (about 0.020742 here).
# imean(income) -> gamma: BETWEEN MINUS WITHIN (about 0.006579 here).
# A heading saying "between" does not make gamma the between slope in this
# parameterisation. The between slope is beta + gamma, about 0.027321.
# Holding mean income fixed, the income coefficient describes within changes.
# H0: gamma = 0 is the Mundlak restriction examined in Lab 2.6.

# B. The within-between specification: equivalent algebra, different labels.
#
# liquor_it = delta + beta_W * (income_it - mean_income_i)
#                   + beta_B * mean_income_i + c_i + e_it
within_between_panelr <- panelr::wbm(
  liquor ~ income,
  data = liquor_panel,
  model = "w-b"
)
summary(within_between_panelr)

# income        -> beta_W: the within slope (about 0.020742 here).
# imean(income) -> beta_B: the between slope (about 0.027321 here).
# Relationships: beta_W = beta; beta_B = beta + gamma.
# Therefore gamma = beta_B - beta_W.
# H0 is now beta_B = beta_W, NOT beta_B = 0. The p-value printed for the
# between coefficient alone tests a zero between slope, not the RE restriction.
# Use the contextual parameterisation to express the restriction as one
# coefficient equal to zero, or test the within-between coefficient contrast.

# Same decomposition, different estimation and inference:
# - Our lm() in Lab 2.5 fits the augmented mean equation using OLS.
# - wbm() fits a mixed model with a household random intercept via lme4.
# - In this balanced, linear, one-regressor example the slope estimates should
#   agree up to numerical precision. Do not extend this claim automatically
#   to random-slope models, nonlinear models or different analysis samples.
# - Default wbm() inference is model-based, with Satterthwaite degrees of
#   freedom. It is NOT our household-clustered CR1 / t(39) convention.
# A random intercept models a particular within-household covariance pattern;
# it does not itself provide arbitrary within-household cluster robustness.
# Consequently, do not expect the default contextual p-value to equal 0.7535.
# Keep the following as the lecture's cluster-robust Mundlak diagnostic:
mundlak_test

# Discussion questions:
# 1. Why does the mean-income coefficient change between the two summaries?
#    Its meaning changes: gamma versus beta + gamma, not a new economic result.
# 2. Does a significant between slope reject the RE restriction?
#    No. The relevant question is whether between and within slopes differ.
# 3. Does adding household means establish a causal income effect?
#    No. Strict exogeneity of the remaining time-varying error is still needed.
# 4. Why might p-values differ even if the slope estimates agree?
#    The covariance assumptions, SE estimators and reference degrees of freedom
#    differ. Compare inference conventions before interpreting the difference.

# Next topics (§15.4.6): common time effects, remaining endogeneity and panel IV,
# and dynamic panels. Household FE alone does not resolve these issues.
