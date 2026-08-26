# IV/2SLS in Practice: Can We Trust the Instrument?
# Standalone live-lab script for SOK-3020 Econometrics
#
# Textbook basis: Principles of Econometrics, 5th ed., Chapter 10,
# Sections 10.3.7--10.4.3.
#
# The script uses base R for the core calculations so students can see exactly
# what each diagnostic is asking. It loads the Mroz data from a local CSV when
# available, from the wooldridge package when installed, or from Rdatasets as a
# final fallback.
#
# IMPORTANT: The numbered LAB STEP sections match the Quarto walkthrough.

# =============================================================================
# SETUP — Load the Mroz data and define helper functions
# =============================================================================

# Clean the environment
rm(list = ls())

# Load the data directly from the POE5 textbook's URL
message("Loading Mroz data from the POE5 URL...")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mroz.rdata"))

# Create log wage safely (handling NA for non-positive wages) and experience squared
mroz$lwage <- ifelse(mroz$wage > 0, log(mroz$wage), NA_real_)
mroz$expersq <- mroz$exper^2       

# Chapter 10 uses the 428 women in the labour force (lfp == 1)
# Subsetting using base R to remain package-light
mroz_work <- subset(mroz, lfp == 1 & is.finite(lwage))
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

# Homoskedastic 2SLS using matrix formulas.
# X contains the structural regressors; Z contains all instruments, including
# the exogenous regressors that instrument themselves.
two_stage_least_squares <- function(y, X, Z) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- nrow(X)
  k <- ncol(X)

  ztz_inv <- solve(crossprod(Z))
  xt_pz_x <- crossprod(X, Z) %*% ztz_inv %*% crossprod(Z, X)
  xt_pz_y <- crossprod(X, Z) %*% ztz_inv %*% crossprod(Z, y)

  beta <- solve(xt_pz_x, xt_pz_y)
  residuals <- as.numeric(y - X %*% beta)
  sigma2 <- sum(residuals^2) / (n - k)
  vcov <- sigma2 * solve(xt_pz_x)
  se <- sqrt(diag(vcov))

  list(
    coefficients = as.numeric(beta),
    se = se,
    vcov = vcov,
    residuals = residuals,
    sigma2 = sigma2,
    n = n,
    k = k
  )
}

# Joint excluded-instrument F test via restricted and unrestricted first stages.
first_stage_f_test <- function(restricted, unrestricted) {
  rss_r <- deviance(restricted)
  rss_u <- deviance(unrestricted)
  q <- length(coef(unrestricted)) - length(coef(restricted))
  df_u <- df.residual(unrestricted)

  F_stat <- ((rss_r - rss_u) / q) / (rss_u / df_u)
  p_value <- pf(F_stat, q, df_u, lower.tail = FALSE)
  partial_r2 <- (rss_r - rss_u) / rss_r

  c(F = F_stat, df1 = q, df2 = df_u, p_value = p_value, partial_R2 = partial_r2)
}

cat("\nMroz labour-force sample size:", nrow(mroz_work), "\n")

# =============================================================================
# LAB STEP 1 — Re-establish the competing OLS and IV estimates
# =============================================================================

# QUESTION: How much does the estimated return to education change when we use
# mother's and father's education as instruments?
#
# PREDICTION STOP: If omitted ability creates upward OLS bias, expect the IV
# estimate to be smaller than OLS.

ols_wage <- lm(lwage ~ educ + exper + expersq, data = mroz_work)

# 1. Define X (Endogenous variable + exogenous controls)
# model.matrix automatically adds the intercept (column of 1s)
X <- model.matrix(~ educ + exper + expersq, data = mroz_work)

# 2. Define Z (Instrument + exogenous controls)
# Replace 'educ' with 'mothereduc', keep the exact same controls
Z <- model.matrix(~ mothereduc + exper + expersq, data = mroz_work)

# 3. Run your custom 2SLS function
iv_wage <- two_stage_least_squares(mroz_work$lwage, X, Z)

names(iv_wage$coefficients) <- colnames(X)
names(iv_wage$se) <- colnames(X)

cat("\nOLS wage equation:\n")
print(round(coef(summary(ols_wage)), 4))

cat("\n2SLS coefficients and homoskedastic standard errors:\n")
iv_table <- cbind(
  Estimate = iv_wage$coefficients,
  Std_Error = iv_wage$se,
  t_value = iv_wage$coefficients / iv_wage$se
)
print(round(iv_table, 4))

# Textbook targets (Example 10.5):
# IV/2SLS: intercept 0.0481, EXPER 0.0442, EXPER^2 -0.0009, EDUC 0.0614.
# Standard errors: approximately 0.4003, 0.0134, 0.0004, 0.0314.
# EDUC t-value is approximately 1.96.

# Optional cross-check with an installed IV package. This block is not required.
if (requireNamespace("ivreg", quietly = TRUE)) {
  iv_pkg <- ivreg::ivreg(
    lwage ~ educ + exper + expersq |
      exper + expersq + mothereduc + fathereduc,
    data = mroz_work
  )
  cat("\nOptional ivreg package cross-check:\n")
  print(round(coef(summary(iv_pkg)), 4))
}

# =============================================================================
# LAB STEP 2 — Are the instruments strong?
# =============================================================================

# QUESTION: Do MOTHEREDUC and FATHEREDUC jointly explain EDUC after controlling
# for EXPER and EXPER^2?
#
# PREDICTION STOP: We want a large excluded-instrument F statistic. The textbook
# uses the familiar F > 10 rule of thumb for the one-endogenous-regressor case.

first_stage_restricted <- lm(educ ~ exper + expersq, data = mroz_work)
first_stage <- lm(
  educ ~ exper + expersq + mothereduc + fathereduc,
  data = mroz_work
)

cat("\nUnrestricted first-stage regression:\n")
print(round(coef(summary(first_stage)), 4))

strength <- first_stage_f_test(first_stage_restricted, first_stage)
cat("\nExcluded-instrument strength diagnostics:\n")
print(round(strength, 4))

# Textbook targets (Example 10.6):
# MOTHEREDUC coefficient 0.1576, t = 4.39.
# FATHEREDUC coefficient 0.1895, t = 5.62.
# Joint first-stage F ≈ 55.40.
# Partial R^2 for both excluded instruments ≈ 0.2076.

# =============================================================================
# LAB STEP 3 — Deliberate failure: strong does not mean valid
# =============================================================================

cat("\nDELIBERATE FAILURE:\n")
cat("'The first-stage F statistic is large, therefore the instruments are valid.'\n")
cat("This is wrong. The first stage is evidence about RELEVANCE only.\n")
cat("Validity additionally requires exclusion/exogeneity: cov(z, e) = 0.\n")

# No additional regression is needed. This is an interpretation checkpoint.

# =============================================================================
# LAB STEP 4 — Do we actually need IV? Regression-based Hausman test
# =============================================================================

# QUESTION: Is EDUC correlated with the structural wage error?
#
# Step 1: use the first-stage residual VHAT.
# Step 2: add VHAT to the original wage equation.
# Test H0: coefficient on VHAT = 0.

mroz_work$vhat <- residuals(first_stage)

hausman_aux <- lm(
  lwage ~ educ + exper + expersq + vhat,
  data = mroz_work
)

cat("\nHausman residual-inclusion regression:\n")
print(round(coef(summary(hausman_aux)), 4))

hausman_row <- coef(summary(hausman_aux))["vhat", ]
cat("\nHausman test of EDUC exogeneity (coefficient on VHAT):\n")
print(round(hausman_row, 4))

# Textbook target (Example 10.7):
# VHAT coefficient 0.0582, SE 0.0348, t = 1.6711, p = 0.0954.
# Interpretation: fail to reject at 5%, reject at 10%; evidence is suggestive,
# not overwhelming.

# =============================================================================
# LAB STEP 5 — Can the surplus restriction survive an overidentification test?
# =============================================================================

# With B = 1 endogenous regressor and L = 2 excluded instruments, there is
# L-B = 1 surplus moment condition. Under homoskedasticity the textbook uses a
# Sargan-style N*R^2 test.

mroz_work$iv_resid <- iv_wage$residuals

overid_aux <- lm(
  iv_resid ~ exper + expersq + mothereduc + fathereduc,
  data = mroz_work
)

overid_r2 <- summary(overid_aux)$r.squared
overid_stat <- nrow(mroz_work) * overid_r2
overid_df <- 1L
overid_p <- pchisq(overid_stat, df = overid_df, lower.tail = FALSE)

cat("\nOveridentification auxiliary regression:\n")
cat("R^2 =", round(overid_r2, 6), "\n")
cat("N * R^2 =", round(overid_stat, 4), "\n")
cat("df =", overid_df, "\n")
cat("p-value =", round(overid_p, 4), "\n")
cat("5% chi-square critical value =", round(qchisq(0.95, df = overid_df), 2), "\n")

# Textbook targets (Example 10.7):
# R^2 = 0.000883, N*R^2 = 0.3779, chi-square(1) 5% critical value = 3.84.
# Therefore fail to reject the surplus restriction.
# IMPORTANT: failure to reject is not proof that both instruments are valid.

# =============================================================================
# LAB STEP 6 — Put the evidence together
# =============================================================================

ols_educ <- coef(summary(ols_wage))["educ", "Estimate"]
iv_educ <- iv_wage$coefficients["educ"]
iv_educ_se <- iv_wage$se["educ"]
hausman_p <- hausman_row["Pr(>|t|)"]

summary_table <- data.frame(
  Evidence = c(
    "OLS EDUC coefficient",
    "IV/2SLS EDUC coefficient",
    "IV/2SLS EDUC standard error",
    "First-stage excluded-instrument F",
    "First-stage partial R-squared",
    "Hausman p-value",
    "Overidentification N*R-squared",
    "Overidentification p-value"
  ),
  Value = c(
    ols_educ,
    iv_educ,
    iv_educ_se,
    strength["F"],
    strength["partial_R2"],
    hausman_p,
    overid_stat,
    overid_p
  )
)

cat("\nEvidence summary:\n")
print(transform(summary_table, Value = round(Value, 4)), row.names = FALSE)

cat("\nInterpretation:\n")
cat("1. The excluded instruments are strongly relevant in the first stage.\n")
cat("2. The Hausman evidence against EDUC exogeneity is suggestive, not decisive.\n")
cat("3. The surplus restriction is not rejected.\n")
cat("4. None of these diagnostics proves the parental-education exclusion assumption.\n")
cat("5. A causal interpretation of the IV coefficient remains conditional on that economic argument.\n")

