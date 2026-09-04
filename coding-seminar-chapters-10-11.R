# =============================================================================
# CODING SEMINAR: ENDOGENEITY, INSTRUMENTAL VARIABLES, AND SIMULTANEOUS EQUATIONS
# Principles of Econometrics, 5th edition — Chapters 10 and 11
# =============================================================================
#
# Seminar purpose
# ----------------
# This script connects two ideas:
#
#   1. Chapter 10: an explanatory variable may be endogenous, so OLS does not
#      identify the structural effect we want.
#   2. Chapter 11: endogeneity may arise naturally because economic variables
#      such as price and quantity are jointly determined in equilibrium.
#
# We use two empirical applications:
#
#   Part A: Returns to education using the MROZ data.
#   Part B: Supply and demand for truffles.
#
# The Fulton Fish Market is deliberately omitted because it was covered in the
# lecture. Seemingly Unrelated Regression (SUR) is also deliberately omitted.
#
# How to use the script
# ---------------------
# Work through the sections in order. At every "PREDICTION STOP", pause before
# running the next command. Discuss what you expect to see and why.
#
# The important questions are not merely "What coefficient did R produce?"
# Instead ask:
#
#   * What is endogenous, and why?
#   * What variation identifies the coefficient?
#   * Is the proposed instrument relevant?
#   * Is its exclusion restriction economically credible?
#   * What conclusion is justified by the evidence?
#
# =============================================================================
# 0. SETUP
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ivreg)

options(digits = 4)

# =============================================================================
# PART A — CHAPTER 10: RETURNS TO EDUCATION AND INSTRUMENTAL VARIABLES
# =============================================================================

# Economic question
# -----------------
# Does another year of education increase a woman's wage, and can parental
# education provide credible exogenous variation in her years of schooling?

# ----------------------------------------------------------------------------
# A1. Load and inspect the MROZ data
# ----------------------------------------------------------------------------

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/mroz.def")
# mroz.def
# 
# taxableinc federaltax hsiblings hfathereduc hmothereduc siblings lfp hours kidsl6 kids618 age educ wage wage76 hhours hage heduc hwage faminc mtr mothereduc fathereduc unemployment largecity exper
# 
# Obs: 753
# 
# variable   		description
# -----------------------------------------------------------------------------
# taxableinc     	Taxable income for household 
# federaltax     	Federal income taxes 
# hsiblings       husband's number of siblings  
# hfathereduc     husband's father's education level
# hmothereduc     husband's mothers's education level
# siblings       	Wife's number of siblings
# lfp        	dummy variable = 1 if woman worked in 1975, else 0
# hours       	Wife's hours of work in 1975
# kidsl6        	Number of children less than 6 years old in household
# kids618       	Number of children between ages 6 and 18 in household
# age         	Wife's age
# educ         	Wife's educational attainment, in years
# wage         	Wife's 1975 average hourly earnings, in 1975 dollars
# wage76       	Wife's wage reported at 1976 interview, for 1976
# hhours       	Husband's hours worked in 1975
# hage         	Husband's age
# heduc         	Husband's educational attainment, in years
# hwage         	Husband's wage, in 1975 dollars
# faminc     	Family income, in 1975 dollars
# mtr        	marginal tax rate facing the wife, includes Soc Sec taxes
# mothereduc      wife's mother's education level
# fathereduc      wife's father's education level
# unemployment   	Unemployment rate in county of residence
# bigcity        	Dummy variable = 1 if live in large city (SMSA), else 0
# exper         	Actual years of wife's previous labor market experience
# 
# 
# THE MROZ DATA FILE IS TAKEN FROM THE 1976 PANEL STUDY OF INCOME
# DYNAMICS, AND IS BASED ON DATA FOR THE PREVIOUS YEAR, 1975.  OF THE 753
# OBSERVATIONS, THE FIRST 428 ARE FOR WOMEN WITH POSITIVE HOURS
# WORKED IN 1975, WHILE THE REMAINING 325 OBSERVATIONS ARE FOR WOMEN
# WHO DID NOT WORK FOR PAY IN 1975.

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mroz.rdata"))

# The wage equation can only be estimated for women with an observed positive
# wage. lfp == 1 denotes labour-force participation.
work <- mroz |>
  filter(lfp == 1, wage > 0)

glimpse(work)

work |>
  summarise(
    observations = n(),
    mean_wage = mean(wage),
    mean_education = mean(educ),
    mean_experience = mean(exper),
    mean_mother_education = mean(mothereduc),
    mean_father_education = mean(fathereduc)
  )

# PREDICTION STOP:
# Before plotting, predict the sign of the unconditional association between
# education and log wages. Does that sign alone have a causal interpretation?

ggplot(work, aes(x = educ, y = log(wage))) +
  geom_jitter(alpha = 0.30, width = 0.15, height = 0) +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
  labs(
    title = "Education and log wages",
    subtitle = "A positive association does not by itself establish causality",
    x = "Years of education",
    y = "Log hourly wage"
  ) +
  theme_minimal()

# ----------------------------------------------------------------------------
# A2. OLS benchmark
# ----------------------------------------------------------------------------

# Experience and experience squared allow the wage profile to rise at a
# decreasing rate over a person's working life.
ols_wage <- lm(
  log(wage) ~ educ + exper + I(exper^2),
  data = work
)

summary(ols_wage)

# In a log-level model, the education coefficient is approximately the
# percentage wage difference associated with one additional year of education.
ols_return <- coef(ols_wage)["educ"]
100 * ols_return
100 * (exp(ols_return) - 1)  # Exact percentage interpretation

# PREDICTION STOP:
# Name at least two reasons why education might be correlated with the wage
# equation error. Examples could include ability, motivation, family background,
# or measurement error. Predict the likely direction of OLS bias—but be prepared
# to explain why the direction is not guaranteed.

# ----------------------------------------------------------------------------
# A3. Candidate instruments and the identification argument
# ----------------------------------------------------------------------------

# Proposed excluded instruments:
#   mothereduc = mother's years of education
#   fathereduc = father's years of education
#
# Relevance question:
# Do these variables help explain the woman's own education?
#
# Exclusion question:
# Conditional on experience, do they affect her wage only through her education?
# This is an economic assumption. It cannot be proved by a first-stage F-test.

work |>
  summarise(
    cor_educ_mother = cor(educ, mothereduc, use = "complete.obs"),
    cor_educ_father = cor(educ, fathereduc, use = "complete.obs")
  )

# Pairwise correlations are descriptive, but instrument relevance in a multiple
# regression is a conditional and joint question. We therefore estimate the
# complete first stage.

# ----------------------------------------------------------------------------
# A4. First stage: do the excluded instruments shift education?
# ----------------------------------------------------------------------------

first_stage_restricted <- lm(
  educ ~ exper + I(exper^2),
  data = work
)

first_stage <- lm(
  educ ~ exper + I(exper^2) + mothereduc + fathereduc,
  data = work
)

summary(first_stage)

# Joint partial F-test for the excluded instruments:
# H0: coefficient(mothereduc) = coefficient(fathereduc) = 0
first_stage_test <- anova(first_stage_restricted, first_stage)
first_stage_test

first_stage_F <- first_stage_test$F[2]
first_stage_p <- first_stage_test$`Pr(>F)`[2]

c(first_stage_F = first_stage_F, first_stage_p = first_stage_p)

# INTERPRETATION:
# A larger partial F-statistic is evidence for relevance. It is not evidence that
# parental education satisfies the exclusion restriction.

# ----------------------------------------------------------------------------
# A5. Manual two-stage calculation: useful for intuition, not inference
# ----------------------------------------------------------------------------

# Stage 1 isolates the part of education predicted by all exogenous variables.
work <- work |>
  mutate(educ_hat = fitted(first_stage))

# PREDICTION STOP:
# Should educ_hat contain the same endogenous variation as observed educ?
# Why should replacing educ with educ_hat alter the estimated wage return?

manual_second_stage <- lm(
  log(wage) ~ educ_hat + exper + I(exper^2),
  data = work
)

summary(manual_second_stage)

# IMPORTANT:
# The coefficient on educ_hat equals the 2SLS coefficient below (up to numerical
# precision), but lm() reports the wrong standard error for 2SLS. Never use this
# manual second-stage standard error for inference.

# ----------------------------------------------------------------------------
# A6. Estimate the wage equation correctly with 2SLS
# ----------------------------------------------------------------------------

# Formula structure:
#   structural regressors | all exogenous variables and excluded instruments
#
# Experience terms appear on both sides because they are included exogenous
# regressors. Parental education appears only after | because it is excluded
# from the structural wage equation.
iv_wage <- ivreg(
  log(wage) ~ educ + exper + I(exper^2) |
    exper + I(exper^2) + mothereduc + fathereduc,
  data = work
)

summary(iv_wage, diagnostics = TRUE)

# Verify the central computational result: manual and dedicated 2SLS give the
# same coefficient, although their reported standard errors differ.
c(
  manual_2sls = coef(manual_second_stage)["educ_hat"],
  ivreg_2sls = coef(iv_wage)["educ"]
)

stopifnot(
  isTRUE(all.equal(
    unname(coef(manual_second_stage)["educ_hat"]),
    unname(coef(iv_wage)["educ"]),
    tolerance = 1e-8
  ))
)

# ----------------------------------------------------------------------------
# A7. Compare OLS and IV as economic estimates
# ----------------------------------------------------------------------------

wage_comparison <- tibble(
  estimator = c("OLS", "2SLS/IV"),
  education_coefficient = c(
    coef(ols_wage)["educ"],
    coef(iv_wage)["educ"]
  ),
  standard_error = c(
    coef(summary(ols_wage))["educ", "Std. Error"],
    coef(summary(iv_wage))["educ", "Std. Error"]
  )
) |>
  mutate(
    approximate_percent = 100 * education_coefficient,
    exact_percent = 100 * (exp(education_coefficient) - 1)
  )

wage_comparison

ggplot(
  wage_comparison,
  aes(x = estimator, y = education_coefficient)
) +
  geom_point(size = 3, colour = "navy") +
  geom_errorbar(
    aes(
      ymin = education_coefficient - 1.96 * standard_error,
      ymax = education_coefficient + 1.96 * standard_error
    ),
    width = 0.12
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Estimated return to education",
    subtitle = "Points are estimates; bars are approximate 95% confidence intervals",
    x = NULL,
    y = "Coefficient in the log-wage equation"
  ) +
  theme_minimal()

# DISCUSSION:
#   1. Why is the IV estimate less precise?
#   2. Does a difference between OLS and IV prove that OLS is biased?
#   3. Which estimate would you report, and what qualification would accompany it?

# ----------------------------------------------------------------------------
# A8. Control-function test for endogeneity
# ----------------------------------------------------------------------------

# Under the null hypothesis that education is exogenous, the unexplained part of
# education from the first stage should not help explain wages. Add the first-
# stage residual to the structural equation and test its coefficient.
work <- work |>
  mutate(educ_residual = residuals(first_stage))

control_function_wage <- lm(
  log(wage) ~ educ + exper + I(exper^2) + educ_residual,
  data = work
)

summary(control_function_wage)

# H0: coefficient(educ_residual) = 0  -> education can be treated as exogenous
# H1: coefficient(educ_residual) != 0 -> evidence that education is endogenous

# Do not interpret failure to reject H0 as proof of exogeneity. The test may have
# limited power, particularly when instruments are weak.

# ----------------------------------------------------------------------------
# A9. Diagnostics from the IV model
# ----------------------------------------------------------------------------

iv_wage_diagnostics <- summary(iv_wage, diagnostics = TRUE)$diagnostics
iv_wage_diagnostics

# The diagnostic table normally includes:
#
#   * Weak-instrument test: addresses relevance.
#   * Wu-Hausman test: compares the exogenous and endogenous specifications.
#   * Sargan test: available because the model is overidentified; evaluates the
#     joint orthogonality restrictions under its assumptions.
#
# A non-rejection in an overidentification test does not prove that the instruments
# are valid. Invalid instruments can fail in similar ways, and tests can have low
# power. Economic reasoning remains essential.

# SEMINAR DECISION — PART A:
# Write two sentences:
#   1. What do the estimates suggest about the return to education?
#   2. What is the most important threat to interpreting the IV result causally?

# =============================================================================
# PART B — CHAPTER 11: SIMULTANEOUS SUPPLY AND DEMAND FOR TRUFFLES
# =============================================================================

# Economic question
# -----------------
# Can we separately estimate the demand and supply responses to truffle price
# when market price and quantity are jointly determined?

# Structural system
# -----------------
# Demand: q = alpha_0 + alpha_p*p + alpha_ps*ps + alpha_di*di + e_d
# Supply: q = beta_0  + beta_p*p  + beta_pf*pf            + e_s
#
# Endogenous variables:
# p = price of premium truffles, $ per ounce
# q = quantity of truffles traded in a market period, in ounces
#
# Exogenous variables:
# ps = price of choice truffles (a substitute), $ per ounce
# di = per capita disposable income, in units of $1000 per month
# pf = hourly rental price of truffle pig, $ per hour

# ----------------------------------------------------------------------------
# B1. Load and inspect the truffle data
# ----------------------------------------------------------------------------

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/truffles.def")

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/truffles.rdata"))

glimpse(truffles)

truffles |>
  summarise(
    observations = n(),
    across(c(p, q, ps, di, pf), list(mean = mean, sd = sd))
  )

ggplot(truffles, aes(x = p, y = q)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
  labs(
    title = "Observed truffle prices and quantities",
    subtitle = "Equilibrium observations do not reveal a structural curve by themselves",
    x = "Premium-truffle price ($ per ounce)",
    y = "Quantity traded (ounces)"
  ) +
  theme_minimal()

# PREDICTION STOP:
# Would a positive price-quantity slope necessarily be a supply curve? Which
# shocks or exogenous variables must move equilibrium for that interpretation?

# ----------------------------------------------------------------------------
# B2. DELIBERATE FAILURE: estimate structural equations by OLS
# ----------------------------------------------------------------------------

ols_demand_truffles <- lm(q ~ p + ps + di, data = truffles)
ols_supply_truffles <- lm(q ~ p + pf, data = truffles)

summary(ols_demand_truffles)
summary(ols_supply_truffles)

# Why this is potentially misleading:
# Price is jointly determined with quantity. Demand and supply shocks affect the
# equilibrium price, so p is generally correlated with each structural error.
# The expected sign or a high R-squared does not repair this simultaneity problem.

# ----------------------------------------------------------------------------
# B3. Identification: which variables shift the other curve?
# ----------------------------------------------------------------------------

identification_table <- tibble(
  equation_to_estimate = c("Demand", "Supply"),
  included_exogenous_controls = c("ps, di", "pf"),
  excluded_instruments = c("pf", "ps, di"),
  economic_shift = c(
    "Factor cost shifts supply and traces demand",
    "Substitute price and income shift demand and trace supply"
  )
)

identification_table

# Exclusion arguments:
#   * pf belongs in supply because it is a production cost, but is excluded from
#     consumer demand.
#   * ps and di belong in demand, but are excluded from production supply.
#
# These restrictions come from economic theory. Do not omit a relevant variable
# merely to satisfy the instrument-counting rule.

# ----------------------------------------------------------------------------
# B4. Reduced form / first stage for the endogenous price
# ----------------------------------------------------------------------------

# Every exogenous variable in the system enters the price reduced form.
price_reduced_form <- lm(p ~ ps + di + pf, data = truffles)
summary(price_reduced_form)

# Equation-specific relevance for demand:
# Is pf relevant for price after controlling for ps and di?
demand_first_stage_restricted <- lm(p ~ ps + di, data = truffles)
demand_relevance_test <- anova(
  demand_first_stage_restricted,
  price_reduced_form
)
demand_relevance_test

# Equation-specific relevance for supply:
# Are ps and di jointly relevant for price after controlling for pf?
supply_first_stage_restricted <- lm(p ~ pf, data = truffles)
supply_relevance_test <- anova(
  supply_first_stage_restricted,
  price_reduced_form
)
supply_relevance_test

# PREDICTION STOP:
# It is possible for the same price reduced form to provide strong identification
# for one structural equation and weak identification for another. Explain why
# before interpreting these two tests.

# The quantity reduced form is useful for describing how equilibrium quantity
# responds to all exogenous market shifters. Its coefficients are equilibrium
# effects, not structural demand or supply coefficients.
quantity_reduced_form <- lm(q ~ ps + di + pf, data = truffles)
summary(quantity_reduced_form)

# ----------------------------------------------------------------------------
# B5. Estimate truffle demand and supply by 2SLS
# ----------------------------------------------------------------------------

# Demand:
#   p is endogenous.
#   ps and di are included exogenous controls.
#   pf is the excluded instrument that shifts supply.
iv_demand_truffles <- ivreg(
  q ~ p + ps + di | ps + di + pf,
  data = truffles
)

# Supply:
#   p is endogenous.
#   pf is an included exogenous production-cost control.
#   ps and di are excluded instruments that shift demand.
iv_supply_truffles <- ivreg(
  q ~ p + pf | pf + ps + di,
  data = truffles
)

summary(iv_demand_truffles, diagnostics = TRUE)
summary(iv_supply_truffles, diagnostics = TRUE)

# INTERPRETATION CHECK:
#   Demand price coefficient: expected to be negative.
#   Substitute-price coefficient: expected to be positive.
#   Income coefficient: expected to be positive for a normal good.
#   Supply price coefficient: expected to be positive.
#   Factor-cost coefficient: expected to be negative.
#
# Do not stop at whether each sign matches theory. Discuss magnitude, uncertainty,
# first-stage strength, and the credibility of the exclusion restrictions.

# ----------------------------------------------------------------------------
# B6. Compare OLS and 2SLS price coefficients
# ----------------------------------------------------------------------------

truffle_price_comparison <- tibble(
  equation = rep(c("Demand", "Supply"), each = 2),
  estimator = rep(c("OLS", "2SLS/IV"), times = 2),
  price_coefficient = c(
    coef(ols_demand_truffles)["p"],
    coef(iv_demand_truffles)["p"],
    coef(ols_supply_truffles)["p"],
    coef(iv_supply_truffles)["p"]
  ),
  standard_error = c(
    coef(summary(ols_demand_truffles))["p", "Std. Error"],
    coef(summary(iv_demand_truffles))["p", "Std. Error"],
    coef(summary(ols_supply_truffles))["p", "Std. Error"],
    coef(summary(iv_supply_truffles))["p", "Std. Error"]
  )
)

truffle_price_comparison

ggplot(
  truffle_price_comparison,
  aes(x = estimator, y = price_coefficient, colour = equation)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  geom_errorbar(
    aes(
      ymin = price_coefficient - 1.96 * standard_error,
      ymax = price_coefficient + 1.96 * standard_error
    ),
    width = 0.12,
    position = position_dodge(width = 0.25)
  ) +
  facet_wrap(~ equation, scales = "free_y") +
  labs(
    title = "Truffle price coefficients: OLS versus 2SLS",
    subtitle = "The estimators use different sources of price variation",
    x = NULL,
    y = "Estimated price coefficient",
    colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ----------------------------------------------------------------------------
# B7. Translate linear coefficients into elasticities at the sample means
# ----------------------------------------------------------------------------

# Because the truffle equations are linear in levels, the price coefficient is
# not itself an elasticity. For y = b*x + ..., elasticity at the means is
# b * mean(x) / mean(y).

truffle_means <- truffles |>
  summarise(across(c(p, q, ps, di, pf), mean))

elasticities_at_means <- tibble(
  equation = c("Demand", "Demand", "Demand", "Supply", "Supply"),
  elasticity = c(
    "Price",
    "Substitute price",
    "Income",
    "Price",
    "Factor cost"
  ),
  value = c(
    coef(iv_demand_truffles)["p"] * truffle_means$p / truffle_means$q,
    coef(iv_demand_truffles)["ps"] * truffle_means$ps / truffle_means$q,
    coef(iv_demand_truffles)["di"] * truffle_means$di / truffle_means$q,
    coef(iv_supply_truffles)["p"] * truffle_means$p / truffle_means$q,
    coef(iv_supply_truffles)["pf"] * truffle_means$pf / truffle_means$q
  )
)

elasticities_at_means

# DISCUSSION:
#   1. Is truffle demand elastic or inelastic at the sample means?
#   2. Are truffles estimated to be normal goods?
#   3. Does the factor-cost elasticity have the economically expected sign?
#   4. Would these elasticities remain constant at different prices and quantities?

# SEMINAR DECISION — PART B:
# Write a short recommendation to a market analyst. State which supply and demand
# responses are supported by the estimates, identify the variation used to
# estimate each equation, and name the most important threat to the conclusion.

# =============================================================================
# FINAL SYNTHESIS
# =============================================================================

# Complete this comparison in words:
#
# MROZ wage model
#   Endogenous regressor:
#   Excluded instruments:
#   Relevance evidence:
#   Main exclusion concern:
#
# Truffle demand
#   Endogenous regressor:
#   Excluded instrument:
#   Economic source of identifying variation:
#   Main exclusion concern:
#
# Truffle supply
#   Endogenous regressor:
#   Excluded instruments:
#   Economic source of identifying variation:
#   Main exclusion concern:
#
# The central lesson:
# A 2SLS command is easy to generate. Credible structural estimation depends on
# explaining why the excluded variables generate relevant variation and why they
# do not enter the structural error through another channel.
