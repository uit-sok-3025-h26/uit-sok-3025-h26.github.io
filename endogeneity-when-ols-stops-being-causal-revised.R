# Endogeneity: When OLS Stops Being Causal
#
# Textbook basis: Principles of Econometrics, 5th ed., Chapter 10,
# Sections 10.1--10.3.6.
#
# The script is deliberately package-light. IV estimation for the simple
# model is implemented directly from the textbook formulas so no IV package is
# required.
#
# IMPORTANT: The numbered LAB STEP sections match the lecture walkthrough.

# =============================================================================
# SETUP — Load the Mroz data, set directory, and define helper functions
# =============================================================================

rm(list = ls())

# 1. Dynamically find the path of this exact script
library(rstudioapi)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 2. Set the working directory to that path
setwd(script_dir)

library(tidyverse)

# This data frame contains data on wages of married and unmarried women.
#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/mroz.def")

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mroz.rdata"))
head(mroz)

# Create log wage
mroz <- 
  mroz |> 
  mutate(lwage = log(wage),
         expersq = exper^2)       

library(AER)

# Chapter 10 uses the 428 women in the labour force.
mroz_work <- mroz |> filter(lfp == 1)
stopifnot(nrow(mroz_work) == 428)

# Simple IV estimator and textbook homoskedastic standard error.
simple_iv <- function(y, x, z) {
  keep <- complete.cases(y, x, z)
  y <- y[keep]
  x <- x[keep]
  z <- z[keep]
  n <- length(y)

  xbar <- mean(x)
  ybar <- mean(y)
  zbar <- mean(z)

  s_zy <- sum((z - zbar) * (y - ybar))
  s_zx <- sum((z - zbar) * (x - xbar))
  s_zz <- sum((z - zbar)^2)

  beta2 <- s_zy / s_zx
  beta1 <- ybar - beta2 * xbar
  resid <- y - beta1 - beta2 * x
  sigma2_iv <- sum(resid^2) / (n - 2)
  var_beta2 <- sigma2_iv * s_zz / (s_zx^2)

  list(
    coefficients = c(intercept = beta1, educ = beta2),
    se_educ = sqrt(var_beta2),
    residuals = resid,
    n = n
  )
}

# Deterministic pseudo-random generator used for the inconsistency simulation.
# All intermediate values remain below 2^53 and are represented exactly by R
# doubles. qnorm() converts the deterministic uniforms to standard normals.
lcg_uniform <- function(n, seed) {
  modulus <- 2^32
  multiplier <- 1664525
  increment <- 1013904223
  state <- seed %% modulus
  out <- numeric(n)

  for (i in seq_len(n)) {
    state <- (multiplier * state + increment) %% modulus
    out[i] <- state / modulus
  }

  out
}

simulate_endogeneity <- function(n, seed = 104729) {
  u <- lcg_uniform(3L * n, seed)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  draws <- matrix(qnorm(u), nrow = 3L, byrow = TRUE)

  ability <- draws[1, ]
  educ_shock <- draws[2, ]
  wage_shock <- draws[3, ]

  # True structural effect of education = 0.08.
  educ <- 12 + 0.8 * ability + educ_shock
  lwage <- 1 + 0.08 * educ + 0.30 * ability + wage_shock

  fit <- lm(lwage ~ educ)
  c(
    n = n,
    beta_hat = unname(coef(fit)["educ"]),
    se = unname(summary(fit)$coefficients["educ", "Std. Error"])
  )
}

figure_dir <- "figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

cat("\nMroz labour-force sample size:", nrow(mroz_work), "\n")

# =============================================================================
# LAB STEP 1 — Estimate the apparently convincing OLS result
# =============================================================================

# QUESTION: How large is the education-wage association after controlling for
# experience and experience squared?
#
# PREDICTION STOP: Expect a positive, precisely estimated EDUC coefficient.

ols_wage <- lm(lwage ~ educ + exper + expersq, data = mroz_work)
print(summary(ols_wage))

anchor_coef <- coef(summary(ols_wage))["educ", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
cat("\nAnchor OLS education result:\n")
print(round(anchor_coef, 4))

# Textbook target (Example 10.1): EDUC about 0.1075, SE about 0.0141.

# =============================================================================
# LAB STEP 2 — Predict the direction of endogeneity bias
# =============================================================================

# No regression is needed here. The reasoning is the econometrics:
#   ability -> education  (+)
#   ability -> wage       (+)
# If ability is omitted, it enters e, giving cov(EDUC, e) > 0.
# Therefore the OLS education coefficient is expected to be biased upward.

cat("\nBias reasoning:\n")
cat("If omitted ability raises both education and wages, cov(EDUC, e) > 0,\n")
cat("so the OLS coefficient is expected to overstate the causal return.\n")

# =============================================================================
# LAB STEP 3 — Deliberate failure: Does more data fix endogeneity?
# =============================================================================

# DGP:
#   EDUC = 12 + 0.8*ABILITY + u
#   LWAGE = 1 + 0.08*EDUC + 0.30*ABILITY + epsilon
# We estimate LWAGE ~ EDUC and deliberately omit ABILITY.
#
# Probability limit:
# beta_OLS -> 0.08 + 0.8*0.30/(0.8^2 + 1) = 0.226341...

true_beta <- 0.08
plim_ols <- true_beta + (0.8 * 0.30) / (0.8^2 + 1)
sample_sizes <- c(200L, 2000L, 20000L, 200000L)

sim_table <- do.call(
  rbind,
  lapply(sample_sizes, simulate_endogeneity, seed = 104729)
)
sim_table <- as.data.frame(sim_table)

cat("\nEndogeneity simulation:\n")
print(transform(sim_table,
  beta_hat = round(beta_hat, 4),
  se = round(se, 4)
))
cat("True structural beta:", true_beta, "\n")
cat("Asymptotic OLS target under omitted ability:", round(plim_ols, 6), "\n")

# Verified deterministic targets (minor last-decimal differences across qnorm
# implementations are possible):
# N=200      beta_hat ~ 0.2291, SE ~ 0.0564
# N=2,000    beta_hat ~ 0.2166, SE ~ 0.0180
# N=20,000   beta_hat ~ 0.2262, SE ~ 0.0057
# N=200,000  beta_hat ~ 0.2271, SE ~ 0.0018
# plim        = 0.226341

png(file.path(figure_dir, "endogeneity-convergence.png"), width = 1280, height = 768, res = 150)
plot(
  sim_table$n, sim_table$beta_hat,
  log = "x", pch = 19,
  xlab = "Sample size (log scale)",
  ylab = "OLS coefficient on education",
  main = "More data make OLS more precise around the wrong target"
)
abline(h = true_beta, lty = 2, lwd = 2)
abline(h = plim_ols, lwd = 2)
legend(
  "right",
  legend = c("OLS estimate", "True causal effect = 0.08", "OLS probability limit ≈ 0.226"),
  pch = c(19, NA, NA),
  lty = c(NA, 2, 1),
  bty = "n"
)
dev.off()

# =============================================================================
# LAB STEP 4 — Does mother's education provide relevant variation?
# =============================================================================

# To introduce IV cleanly, follow Example 10.2 and temporarily simplify to:
#   LWAGE = beta1 + beta2*EDUC + e
# Instrument: MOTHEREDUC.

ols_simple <- lm(lwage ~ educ, data = mroz_work)
first_stage <- lm(educ ~ mothereduc, data = mroz_work)

cat("\nSimple OLS wage equation:\n")
print(round(coef(summary(ols_simple)), 4))

cat("\nRaw correlation between MOTHEREDUC and EDUC:\n")
print(cor(mroz_work$mothereduc, mroz_work$educ))

cat("\nFirst-stage regression EDUC ~ MOTHEREDUC:\n")
print(round(coef(summary(first_stage)), 4))

# Textbook targets (Examples 10.2 and 10.3):
# Simple OLS EDUC coefficient ~ 0.1086, SE ~ 0.0144.
# First stage: intercept ~ 10.1145, MOTHEREDUC ~ 0.2674, t ~ 8.66.
# The chapter reports a raw EDUC-MOTHEREDUC correlation of about 0.387.

# =============================================================================
# LAB STEP 5 — Replace OLS variation with instrumental-variable variation
# =============================================================================

iv_simple <- with(
  mroz_work,
  simple_iv(y = lwage, x = educ, z = mothereduc)
)

cat("\nSimple IV estimate using MOTHEREDUC as instrument:\n")
print(round(iv_simple$coefficients, 4))
cat("Correct textbook homoskedastic IV SE for EDUC:", round(iv_simple$se_educ, 4), "\n")

# Textbook target (Example 10.2):
# intercept ~ 0.7022, EDUC ~ 0.0385, correct SE ~ 0.0382.

comparison <- data.frame(
  estimator = c("OLS", "IV"),
  educ_coefficient = c(unname(coef(ols_simple)["educ"]), iv_simple$coefficients["educ"]),
  se = c(coef(summary(ols_simple))["educ", "Std. Error"], iv_simple$se_educ)
)

cat("\nOLS versus IV:\n")
print(transform(comparison,
  educ_coefficient = round(educ_coefficient, 4),
  se = round(se, 4)
), row.names = FALSE)

# =============================================================================
# LAB STEP 6 — Deliberate failure: manual 2SLS gives the wrong usual OLS SE
# =============================================================================

mroz_work$educ_hat <- fitted(first_stage)
manual_second_stage <- lm(lwage ~ educ_hat, data = mroz_work)
manual_table <- coef(summary(manual_second_stage))

cat("\nManual second-stage regression LWAGE ~ EDUC_HAT:\n")
print(round(manual_table, 4))

cat("\nCoefficient check:\n")
cat("Manual 2SLS slope:", round(unname(coef(manual_second_stage)["educ_hat"]), 4), "\n")
cat("Direct IV slope:   ", round(unname(iv_simple$coefficients["educ"]), 4), "\n")
cat("Manual OLS SE (WRONG for IV inference):",
    round(manual_table["educ_hat", "Std. Error"], 4), "\n")
cat("Correct IV SE:                         ",
    round(iv_simple$se_educ, 4), "\n")

# Textbook target (Example 10.3): manual second-stage coefficient ~ 0.0385,
# manual ordinary OLS SE ~ 0.0396 (incorrect), correct IV SE ~ 0.0382.

# =============================================================================
# OPTIONAL — Compact results table for lecture preparation
# =============================================================================

results_summary <- data.frame(
  quantity = c(
    "Multiple OLS EDUC coefficient",
    "Multiple OLS EDUC SE",
    "Simple OLS EDUC coefficient",
    "Simple OLS EDUC SE",
    "First-stage MOTHEREDUC coefficient",
    "First-stage MOTHEREDUC t",
    "Simple IV EDUC coefficient",
    "Simple IV EDUC SE",
    "Manual second-stage usual OLS SE (wrong)"
  ),
  value = c(
    coef(ols_wage)["educ"],
    coef(summary(ols_wage))["educ", "Std. Error"],
    coef(ols_simple)["educ"],
    coef(summary(ols_simple))["educ", "Std. Error"],
    coef(first_stage)["mothereduc"],
    coef(summary(first_stage))["mothereduc", "t value"],
    iv_simple$coefficients["educ"],
    iv_simple$se_educ,
    manual_table["educ_hat", "Std. Error"]
  )
)

cat("\nCompact results summary:\n")
print(transform(results_summary, value = round(value, 4)), row.names = FALSE)

# =============================================================================
# END OF LIVE LAB
# =============================================================================

cat("\nTake-home idea:\n")
cat("Precision is not identification. A very precise OLS coefficient can converge\n")
cat("to the wrong causal value when the regressor is endogenous. IV replaces the\n")
cat("contaminated variation in x with variation supplied by a credible instrument.\n")

