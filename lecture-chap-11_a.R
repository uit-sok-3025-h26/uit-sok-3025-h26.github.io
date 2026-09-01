# Lecture Chapter 11A: Simultaneous equations, OLS failure, and identification
#
# This script is deliberately self-contained and uses only base R.
# It matches the numbered live-lab steps in lecture-11a.qmd.

options(digits = 4)

# -----------------------------------------------------------------------------
# 0. Construct equilibrium data
# -----------------------------------------------------------------------------

# The deterministic shocks make every run identical. They behave like irregular
# market shocks without requiring a particular random-number generator.
n <- 240
i <- seq_len(n)

income <- as.numeric(scale(sin(0.37 * i) + cos(0.11 * i)))
demand_shock <- 2 * sin(1.17 * i)
supply_shock <- 2 * cos(0.83 * i)

# Structural demand: Qd = 100 - 2P + 8X + ed
alpha_0 <- 100
alpha_1 <- -2
alpha_2 <- 8

# Structural supply: Qs = 20 + 1.5P + es
beta_0 <- 20
beta_1 <- 1.5

# Solve Qd = Qs for the jointly determined equilibrium price and quantity.
price <- (alpha_0 - beta_0 + alpha_2 * income +
            demand_shock - supply_shock) / (beta_1 - alpha_1)
quantity <- beta_0 + beta_1 * price + supply_shock

market <- data.frame(
  market = i,
  income,
  price,
  quantity,
  demand_shock,
  supply_shock
)

# Confirm that the constructed observations satisfy both structural equations.
quantity_demanded <- alpha_0 + alpha_1 * price +
  alpha_2 * income + demand_shock
stopifnot(max(abs(quantity - quantity_demanded)) < 1e-10)

# -----------------------------------------------------------------------------
# 1. What does the observed price-quantity scatterplot appear to show?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Will the observed relationship be positive or negative? Does its sign alone
# identify a supply curve?

plot(
  quantity ~ price,
  data = market,
  pch = 19,
  col = adjustcolor("steelblue", alpha.f = 0.55),
  xlab = "Equilibrium price",
  ylab = "Equilibrium quantity",
  main = "Observed equilibrium outcomes"
)

naive_supply <- lm(quantity ~ price, data = market)
abline(naive_supply, col = "firebrick", lwd = 2)

# -----------------------------------------------------------------------------
# 2. DELIBERATE FAILURE: treat equilibrium price as exogenous
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# The true supply slope is 1.5. Given that a positive supply shock lowers the
# equilibrium price, should OLS be above or below 1.5?

summary(naive_supply)
coef(naive_supply)
summary(naive_supply)$r.squared
cor(market$price, market$supply_shock)

# Verified targets used in the lecture:
stopifnot(abs(coef(naive_supply)[2] - 1.388705) < 1e-5)
stopifnot(abs(summary(naive_supply)$r.squared - 0.847823) < 1e-5)
stopifnot(abs(cor(market$price, market$supply_shock) + 0.185869) < 1e-5)

# A high R-squared and the expected positive sign do not make price exogenous.

# -----------------------------------------------------------------------------
# 3. What do the reduced forms tell us?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Income shifts demand outward. What should happen to equilibrium price and
# equilibrium quantity after the market adjusts?

price_reduced_form <- lm(price ~ income, data = market)
quantity_reduced_form <- lm(quantity ~ income, data = market)

summary(price_reduced_form)
summary(quantity_reduced_form)

coef(price_reduced_form)
coef(quantity_reduced_form)

# Theoretical equilibrium effect of income on price.
alpha_2 / (beta_1 - alpha_1)

stopifnot(abs(coef(price_reduced_form)[2] - 2.293196) < 1e-5)
stopifnot(abs(coef(quantity_reduced_form)[2] - 3.417255) < 1e-5)

# These are equilibrium effects, not structural supply or demand coefficients.

# -----------------------------------------------------------------------------
# 4. Which curve do income shifts trace?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# If income is excluded from supply, should income-group equilibrium means line
# up with supply or demand?

market$income_group <- cut(
  market$income,
  breaks = quantile(market$income, probs = seq(0, 1, 0.25)),
  include.lowest = TRUE
)

group_means <- aggregate(
  cbind(price, quantity) ~ income_group,
  data = market,
  FUN = mean
)

plot(
  quantity ~ price,
  data = market,
  pch = 19,
  col = adjustcolor("grey40", alpha.f = 0.25),
  xlab = "Equilibrium price",
  ylab = "Equilibrium quantity",
  main = "Demand shifts trace the supply curve"
)
points(
  quantity ~ price,
  data = group_means,
  pch = 19,
  col = "firebrick",
  cex = 1.7
)
abline(a = beta_0, b = beta_1, col = "navy", lwd = 2)
legend(
  "topleft",
  legend = c("Market observations", "Income-group means", "True supply"),
  col = c("grey40", "firebrick", "navy"),
  pch = c(19, 19, NA),
  lty = c(NA, NA, 1),
  bty = "n"
)

# -----------------------------------------------------------------------------
# 5. CORRECTION: isolate income-induced price variation
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Will the coefficient based on predicted price be closer to the true supply
# slope of 1.5 than the naive OLS estimate?

market$predicted_price <- fitted(price_reduced_form)
manual_second_stage <- lm(quantity ~ predicted_price, data = market)
coef(manual_second_stage)

stopifnot(abs(coef(manual_second_stage)[2] - 1.490172) < 1e-5)

# IMPORTANT: the coefficient illustrates the IV logic, but the ordinary OLS
# standard errors from this manual second stage are not valid 2SLS standard
# errors. Use dedicated IV software for inference.

# -----------------------------------------------------------------------------
# 6. DELIBERATE FAILURE: add income and call the result demand
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Does controlling for the demand shifter make observed price exogenous in the
# demand equation?

naive_demand <- lm(quantity ~ price + income, data = market)
summary(naive_demand)
coef(naive_demand)

stopifnot(abs(coef(naive_demand)[2] + 0.244616) < 1e-5)

# The true demand price coefficient is -2, not approximately -0.245. Demand is
# unidentified because the system contains no exogenous supply shifter that is
# excluded from demand.

# -----------------------------------------------------------------------------
# Compact results table used by the lecture page
# -----------------------------------------------------------------------------

results <- data.frame(
  quantity = c(
    "Naive OLS supply slope",
    "True supply slope",
    "Income coefficient in price reduced form",
    "Manual IV/2SLS coefficient",
    "Naive OLS demand price coefficient",
    "True demand price coefficient"
  ),
  estimate = c(
    coef(naive_supply)[2],
    beta_1,
    coef(price_reduced_form)[2],
    coef(manual_second_stage)[2],
    coef(naive_demand)[2],
    alpha_1
  )
)

print(results, row.names = FALSE)

