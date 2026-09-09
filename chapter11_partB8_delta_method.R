# ----------------------------------------------------------------------------
# B8. Delta-method confidence intervals for the elasticities
# ----------------------------------------------------------------------------

# Section B7 reported elasticity point estimates and stopped there. A point
# estimate without an interval invites the reader to treat "elastic" or "normal
# good" as an established fact rather than an estimate. The elasticity is a
# nonlinear function of estimated coefficients, so its standard error must be
# derived rather than read off the regression output.

# The delta method
# ----------------
# If b is asymptotically normal with covariance V, and g(.) is continuously
# differentiable, then
#
#   g(b) ~approx normal( g(beta) , G' V G ),  where G = d g(beta) / d beta.
#
# In words: linearise the nonlinear function around the true parameter and
# propagate the coefficient covariance through that linear approximation. The
# gradient must be evaluated at the estimates, and V must be the covariance
# matrix of the estimator actually used (here 2SLS, not OLS).

# Requires the car package for the deltaMethod() cross-check in B8.4.
# install.packages("car")

# ----------------------------------------------------------------------------
# B8.1 Two different targets, not two different formulas
# ----------------------------------------------------------------------------

# Elasticity at the means can be written in two ways that look identical but
# describe different estimands:
#
#   (i)  FIXED DENOMINATOR
#        eps = b_x * xbar / qbar, treating xbar and qbar as fixed numbers that
#        merely say where we evaluate. The function is then LINEAR in b, so the
#        delta method collapses to rescaling: se(eps) = |xbar/qbar| * se(b_x).
#
#   (ii) FITTED DENOMINATOR
#        eps = b_x * xbar / qhat, where qhat = xbar' b is the quantity the
#        estimated equation predicts at mean regressors. Now b appears in the
#        denominator too, the function is genuinely nonlinear, and every
#        coefficient in the equation contributes to the standard error.
#
# Numerical curiosity worth flagging to students: because the intercept is its
# own instrument, the 2SLS residuals sum to zero, so qhat = qbar exactly. The
# two approaches therefore give the SAME point estimate and DIFFERENT standard
# errors. Nothing has gone wrong; they answer different questions about what is
# being held fixed.

# ----------------------------------------------------------------------------
# B8.2 A hand-coded delta-method function
# ----------------------------------------------------------------------------

# For eps = b_j * xbar_j / (xbar' b) the gradient has a compact closed form:
#
#   d eps / d b = (xbar_j / qhat) * e_j  -  (eps / qhat) * xbar
#
# where e_j is the selector vector for the coefficient of interest. The first
# term is the direct numerator effect; the second is the feedback through the
# fitted denominator. Setting the second term to zero returns case (i).

elasticity_delta <- function(model,
                             regressor,
                             data,
                             level = 0.95,
                             denominator = c("fitted", "sample_mean")) {

  denominator <- match.arg(denominator)

  b <- coef(model)
  V <- vcov(model)

  regressor_names <- names(b)[names(b) != "(Intercept)"]
  xbar <- c(1, colMeans(data[, regressor_names, drop = FALSE]))
  names(xbar) <- names(b)

  qhat <- sum(b * xbar)    # fitted quantity at mean regressors
  qbar <- mean(data$q)     # observed mean quantity

  grad <- rep(0, length(b))
  names(grad) <- names(b)

  if (denominator == "sample_mean") {
    scale_factor <- xbar[[regressor]] / qbar
    estimate <- scale_factor * b[[regressor]]
    grad[regressor] <- scale_factor
  } else {
    estimate <- b[[regressor]] * xbar[[regressor]] / qhat
    grad <- -(estimate / qhat) * xbar
    grad[regressor] <- grad[[regressor]] + xbar[[regressor]] / qhat
  }

  std_error <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  critical_value <- qt(1 - (1 - level) / 2, df = df.residual(model))

  tibble(
    regressor = regressor,
    denominator = denominator,
    elasticity = estimate,
    std_error = std_error,
    conf_low = estimate - critical_value * std_error,
    conf_high = estimate + critical_value * std_error
  )
}

# PREDICTION STOP:
# Before running the next block: which of the two denominators will give the
# wider interval for the demand price elasticity, and why? Is the answer
# guaranteed in general, or does it depend on the covariances between the
# coefficients?

# ----------------------------------------------------------------------------
# B8.3 Intervals for all five elasticities
# ----------------------------------------------------------------------------

elasticity_specification <- tibble(
  equation = c("Demand", "Demand", "Demand", "Supply", "Supply"),
  label = c("Price", "Substitute price", "Income", "Price", "Factor cost"),
  regressor = c("p", "ps", "di", "p", "pf")
)

elasticity_intervals <- elasticity_specification |>
  rowwise() |>
  reframe(
    equation = equation,
    label = label,
    bind_rows(
      elasticity_delta(
        model = if (equation == "Demand") iv_demand_truffles else iv_supply_truffles,
        regressor = regressor,
        data = truffles,
        denominator = "sample_mean"
      ),
      elasticity_delta(
        model = if (equation == "Demand") iv_demand_truffles else iv_supply_truffles,
        regressor = regressor,
        data = truffles,
        denominator = "fitted"
      )
    )
  )

elasticity_intervals |>
  mutate(across(c(elasticity, std_error, conf_low, conf_high), \(x) round(x, 3))) |>
  print(n = Inf)

# Confirm that the B7 point estimates are reproduced exactly.
elasticity_intervals |>
  filter(denominator == "fitted") |>
  select(equation, label, elasticity)

# ----------------------------------------------------------------------------
# B8.4 Cross-check against car::deltaMethod()
# ----------------------------------------------------------------------------

# Never trust a hand-coded gradient without a second opinion. car::deltaMethod()
# differentiates the expression symbolically, so agreement to several decimals
# is evidence that the algebra in B8.2 is right.

# car::deltaMethod() has no ivreg method, so parameterNames is ignored. Rename
# the coefficients and the covariance matrix explicitly instead.

truffle_constants <- c(
  pbar  = truffle_means$p,
  psbar = truffle_means$ps,
  dibar = truffle_means$di,
  pfbar = truffle_means$pf,
  qbar  = mean(truffles$q)
)

rename_parameters <- function(model, new_names) {
  b <- coef(model)
  V <- vcov(model)
  names(b) <- new_names
  dimnames(V) <- list(new_names, new_names)
  list(coefficients = b, vcov = V)
}

demand_pars <- rename_parameters(iv_demand_truffles, c("b0", "bp", "bps", "bdi"))
supply_pars <- rename_parameters(iv_supply_truffles, c("b0", "bp", "bpf"))

# Demand price elasticity, fitted denominator.
car::deltaMethod(
  demand_pars$coefficients,
  "bp * pbar / (b0 + bp * pbar + bps * psbar + bdi * dibar)",
  vcov. = demand_pars$vcov,
  constants = truffle_constants
)

# Demand price elasticity, fixed denominator.
car::deltaMethod(
  demand_pars$coefficients,
  "bp * pbar / qbar",
  vcov. = demand_pars$vcov,
  constants = truffle_constants
)

# Supply price elasticity, fitted denominator.
car::deltaMethod(
  supply_pars$coefficients,
  "bp * pbar / (b0 + bp * pbar + bpf * pfbar)",
  vcov. = supply_pars$vcov,
  constants = truffle_constants
)

# ----------------------------------------------------------------------------
# B8.5 Visualise the intervals
# ----------------------------------------------------------------------------

ggplot(
  elasticity_intervals,
  aes(x = label, y = elasticity, colour = denominator)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -1, linetype = "dotted", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.35)) +
  geom_errorbar(
    aes(ymin = conf_low, ymax = conf_high),
    width = 0.15,
    position = position_dodge(width = 0.35)
  ) +
  coord_flip() +
  facet_wrap(~ equation, scales = "free_y") +
  labs(
    title = "Delta-method 95% intervals for the truffle elasticities",
    subtitle = "Dashed line: zero. Dotted line: unit elasticity.",
    x = NULL,
    y = "Elasticity at the sample means",
    colour = "Denominator"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ----------------------------------------------------------------------------
# B8.6 How good is the linear approximation? A bootstrap comparison
# ----------------------------------------------------------------------------

# The delta method is a first-order approximation justified asymptotically. With
# n = 30 and a ratio of estimates, the sampling distribution can be skewed and
# the symmetric interval can be a poor summary. A pairs bootstrap re-estimates
# the whole 2SLS system on resampled markets and needs no normality assumption.

set.seed(3020)

bootstrap_replicates <- 2000

bootstrap_demand_price_elasticity <- replicate(bootstrap_replicates, {
  resampled <- truffles[sample(nrow(truffles), replace = TRUE), ]

  fit <- try(
    ivreg(q ~ p + ps + di | ps + di + pf, data = resampled),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NA_real_)

  b <- coef(fit)
  xbar <- c(1, colMeans(resampled[, c("p", "ps", "di")]))
  unname(b["p"] * xbar[2] / sum(b * xbar))
})

bootstrap_summary <- tibble(
  method = c("Delta method (fitted)", "Bootstrap percentile"),
  conf_low = c(
    elasticity_intervals |>
      filter(equation == "Demand", label == "Price", denominator == "fitted") |>
      pull(conf_low),
    quantile(bootstrap_demand_price_elasticity, 0.025, na.rm = TRUE)
  ),
  conf_high = c(
    elasticity_intervals |>
      filter(equation == "Demand", label == "Price", denominator == "fitted") |>
      pull(conf_high),
    quantile(bootstrap_demand_price_elasticity, 0.975, na.rm = TRUE)
  )
)

bootstrap_summary

ggplot(
  tibble(elasticity = bootstrap_demand_price_elasticity),
  aes(x = elasticity)
) +
  geom_histogram(bins = 40, alpha = 0.8) +
  geom_vline(
    data = bootstrap_summary,
    aes(xintercept = conf_low, colour = method),
    linetype = "dashed"
  ) +
  geom_vline(
    data = bootstrap_summary,
    aes(xintercept = conf_high, colour = method),
    linetype = "dashed"
  ) +
  labs(
    title = "Bootstrap distribution of the demand price elasticity",
    subtitle = "Compare the symmetric delta-method limits with the percentile limits",
    x = "Elasticity at the means",
    y = "Bootstrap replications",
    colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Better looking plot

bootstrap_window <- c(-10, 1)

share_outside <- mean(
  bootstrap_demand_price_elasticity < bootstrap_window[1] |
    bootstrap_demand_price_elasticity > bootstrap_window[2],
  na.rm = TRUE
)

share_outside

ggplot(
  tibble(elasticity = bootstrap_demand_price_elasticity),
  aes(x = elasticity)
) +
  geom_histogram(binwidth = 0.25, alpha = 0.8) +
  geom_vline(
    data = bootstrap_summary,
    aes(xintercept = conf_low, colour = method),
    linetype = "dashed"
  ) +
  geom_vline(
    data = bootstrap_summary,
    aes(xintercept = conf_high, colour = method),
    linetype = "dashed"
  ) +
  coord_cartesian(xlim = bootstrap_window) +
  labs(
    title = "Bootstrap distribution of the demand price elasticity",
    subtitle = sprintf(
      "View truncated to [%g, %g]; %.1f%% of replications fall outside",
      bootstrap_window[1], bootstrap_window[2], 100 * share_outside
    ),
    x = "Elasticity at the means",
    y = "Bootstrap replications",
    colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# ----------------------------------------------------------------------------
# DISCUSSION — PART B8
# ----------------------------------------------------------------------------
# 1. Does the demand price elasticity interval lie entirely below -1? If it
#    straddles -1, the sample cannot settle whether demand is elastic, even
#    though the point estimate takes a side.
# 2. Does the income elasticity interval exclude zero? State what the interval
#    licenses you to say about truffles being a normal good, and what it does not.
# 3. The delta method inherits every assumption behind V. If the instruments are
#    weak, the 2SLS covariance matrix is unreliable and so is any interval built
#    on it. Revisit the first-stage tests from B4 before quoting these numbers.
# 4. Why is the delta-method interval symmetric by construction, and why is a
#    ratio of estimates a case where symmetry is a substantive assumption rather
#    than a harmless convenience?
# 5. Heteroskedasticity: replace vcov(model) with a robust covariance matrix
#    (sandwich::vcovHC) inside elasticity_delta() and rerun. Which conclusions
#    survive?

# SEMINAR DECISION — PART B8:
# Rewrite the recommendation from Part B so that every elasticity claim is stated
# as an interval rather than a number. Mark each claim the data support, each
# claim the data cannot separate from its alternative, and state which single
# assumption the intervals would collapse without.


