# =============================================================================
# panel_data_01.R
# Panel Data and Fixed Effects · Comparing households with themselves
#
# Accompanies the lecture panel_data_01.
# Hill, Griffiths & Lim, Principles of Econometrics (5th ed.), Chapter 15,
# §§15.1–15.2.4, with an introduction to §§15.3.1–15.3.2.
#
# Run one numbered section at a time. Before you run a lab, write down your
# prediction for the question stated in the lecture.
#
# Blocks marked [EXTRA] go beyond the lecture page. They are optional in
# class and intended for self-study.
# =============================================================================

rm(list = ls())  # clear workspace
getwd()

# Packages --------------------------------------------------------------------

library(tidyverse)
library(plm)       # panel estimators (note: masks dplyr::between, lag, lead)
library(sandwich)  # cluster-robust covariance matrices
library(lmtest)    # coeftest() and coefci()


# Data ------------------------------------------------------------------------
# Definition file:
# https://www.principlesofeconometrics.com/poe5/data/def/liquor5.def
#   hh      household identifier
#   year    period, coded 1–3 (not calendar years)
#   liquor  liquor expenditure, thousands of dollars
#   income  income, thousands of dollars

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/liquor5.rdata"))

panel <- as_tibble(liquor5)
panel

# Is the panel balanced? Every household should appear three times.
count(panel, hh) |> count(n, name = "households")

N  <- n_distinct(panel$hh)    # households (= clusters, G)
TT <- n_distinct(panel$year)  # periods (TT, because T is shorthand for TRUE)
NT <- nrow(panel)             # observations
c(N = N, T = TT, NT = NT)


# 0. Retrieval and connection -------------------------------------------------
# No code. Suppose a regression of liquor on income omits stable household
# preferences. When would that omission bias the income coefficient?
# We return to this question in Lab 1.4 (F-test) and in Section 4 (simulation).


# 1. The hook -----------------------------------------------------------------

pooled <- lm(liquor ~ income, data = panel)
coef(pooled)

# Units: both variables are in $1,000, so the slope is $1,000 of liquor per
# $1,000 of income. In dollars per $1,000 of income:
1000 * coef(pooled)["income"]

# Prediction: will removing between-household differences increase or
# decrease the income coefficient? Write down an economic mechanism.


# 3. The econometric lab ------------------------------------------------------

# Lab 1.1 · Where is the income variation? ------------------------------------

panel <- panel |>
  group_by(hh) |>
  mutate(
    income_mean   = mean(income),
    liquor_mean   = mean(liquor),
    income_within = income - income_mean,
    liquor_within = liquor - liquor_mean
  ) |>
  ungroup()

# Sums of squares in the handout's notation: T = W + B
Txx <- sum((panel$income - mean(panel$income))^2)
Wxx <- sum(panel$income_within^2)
Bxx <- sum((panel$income_mean - mean(panel$income))^2)
c(Txx = Txx, Wxx = Wxx, Bxx = Bxx, W_plus_B = Wxx + Bxx)

within_share <- Wxx / Txx
within_share   # about 3.78 % of the income variation is within households

# For comparison: liquor expenditure varies much more within households
sum(panel$liquor_within^2) / sum((panel$liquor - mean(panel$liquor))^2)

# [EXTRA] One household up close (POE5 Exercise 15.16a)
panel |>
  filter(hh == 1) |>
  select(hh, year, income, income_mean, income_within)

sum(panel$income_within[panel$hh == 1])   # deviations always sum to zero

# Figure: share of income variation within and between households
# (saved as figures/variation_share.png, used in the lecture)
fig_variation_share <- tibble(
  source = factor(c("Within households", "Between households"),
                  levels = c("Within households", "Between households")),
  share  = c(within_share, 1 - within_share)
) |>
  ggplot(aes(x = source, y = share)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::percent(share, accuracy = 0.01)), vjust = -0.4) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(x = NULL, y = "Share of total income sum of squares",
       title = "Where is the income variation?") +
  theme_minimal(base_size = 12)

fig_variation_share

# Figure: levels with pooled fit (left) and within deviations with FE fit (right).
# Orange points are household means; they all sit at the origin after demeaning,
# so they are drawn in the left panel only.
# To trace individual households, set highlight_hh to a few household numbers.
view_levels <- c("Levels: observations, household means, pooled fit",
                 "Deviations from household means: within fit")

plot_obs <- bind_rows(
  panel |> transmute(hh, x = income, y = liquor, view = view_levels[1]),
  panel |> transmute(hh, x = income_within, y = liquor_within, view = view_levels[2])
) |>
  mutate(view = factor(view, levels = view_levels))

plot_means <- panel |>
  distinct(hh, x = income_mean, y = liquor_mean) |>
  mutate(view = factor(view_levels[1], levels = view_levels))

highlight_hh <- integer(0)   # e.g. c(1, 12, 33)

plot_obs <- plot_obs |> mutate(traced = hh %in% highlight_hh)

fig_within_between <- ggplot(plot_obs, aes(x = x, y = y)) +
  geom_point(aes(colour = traced), size = 1.2, alpha = 0.7) +
  geom_path(data = filter(plot_obs, traced), aes(group = hh),
            colour = "firebrick", linewidth = 0.5) +
  geom_point(data = plot_means, shape = 21, fill = "orange", size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = "firebrick"),
                      guide = "none") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  facet_wrap(~ view, scales = "free") +
  labs(x = "Income ($1,000): levels or deviations",
       y = "Liquor expenditure ($1,000): levels or deviations",
       caption = "Axes differ: demeaning changes the scale of the data.") +
  theme_minimal(base_size = 12)

fig_within_between

# Save both figures at the paths the lecture expects. The lecture links to
# ../../figures/, i.e. two levels above the lecture file; adjust figure_dir if
# your script sits elsewhere in the repository.

# Install if you don't have it: install.packages("this.path")
library(this.path)

# Sets directory to the folder containing the script
setwd(this.dir())
getwd()

figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(figure_dir, "variation_share.png"), fig_variation_share,
       width = 6, height = 4, dpi = 300, bg = "white")

ggsave(file.path(figure_dir, "within_between.png"), fig_within_between,
       width = 9, height = 4.5, dpi = 300, bg = "white")


# Lab 1.2 · Which comparison does pooling resemble? ---------------------------

household_means <- panel |>
  distinct(hh, income_mean, liquor_mean) |>
  arrange(hh)

between <- lm(liquor_mean ~ income_mean, data = household_means)

coef(pooled)    # slope 0.02707
coef(between)   # slope 0.02732

# Would adding more households fix the interpretation? More households
# improve precision; they do not change which comparison the slope makes.


# Lab 1.3 · What changes when a household is compared with itself? ------------

within_manual <- lm(liquor_within ~ 0 + income_within, data = panel)

Wxy <- sum(panel$income_within * panel$liquor_within)
c(Wxx = Wxx, Wxy = Wxy, b_FE = Wxy / Wxx)
coef(within_manual)

1000 * coef(within_manual)   # about $20.74 per $1,000 of income

# Note: the default OLS standard errors from within_manual are NOT the final
# panel inference (see Lab 1.4 and Lab 1.5).

# [EXTRA] The pooled slope is an exact weighted average of the within and
# between slopes, with weights given by the variation shares from Lab 1.1:
#   b_pooled = (Wxx/Txx) * b_within + (Bxx/Txx) * b_between
# (exact for a balanced panel). This is why pooling resembles the between
# regression here: 96 % of the weight is on between-household comparisons.
b_pooled  <- coef(pooled)[["income"]]
b_between <- coef(between)[["income_mean"]]
b_within  <- coef(within_manual)[["income_within"]]

c(weighted_average = within_share * b_within + (1 - within_share) * b_between,
  pooled           = b_pooled)

# [EXTRA] FE as a weighted average of household-specific slopes (the handout's
# starting point). Each household's own OLS slope is b_i = Wxy_i / Wxx_i, and
#   b_FE = sum_i (Wxx_i / Wxx) * b_i
# Households with more income variation over time get more weight.
household_slopes <- panel |>
  group_by(hh) |>
  summarise(
    Wxx_i = sum(income_within^2),
    Wxy_i = sum(income_within * liquor_within),
    b_i   = Wxy_i / Wxx_i
  ) |>
  mutate(weight = Wxx_i / sum(Wxx_i))

summary(household_slopes$b_i)   # each slope uses only three observations

c(weighted_mean_b_i   = sum(household_slopes$weight * household_slopes$b_i),
  b_FE                = b_within,
  unweighted_mean_b_i = mean(household_slopes$b_i))

ggplot(household_slopes, aes(x = b_i, y = weight)) +
  geom_point() +
  geom_vline(xintercept = b_within, linetype = "dashed") +
  labs(x = "Household-specific OLS slope (3 observations each)",
       y = "Weight in the FE estimate (Wxx_i / Wxx)",
       title = "The FE slope is a variance-weighted average of household slopes") +
  theme_minimal()

# [EXTRA] First differences versus within (POE5 §15.2.1–15.2.2, Exercise 15.17a)
# dplyr::lag is written in full because plm masks lag().
panel_fd <- panel |>
  arrange(hh, year) |>
  group_by(hh) |>
  mutate(
    d_liquor = liquor - dplyr::lag(liquor),
    d_income = income - dplyr::lag(income)
  ) |>
  ungroup()

first_diff <- lm(d_liquor ~ 0 + d_income, data = panel_fd)  # year 1 rows dropped (NA)
c(FD = coef(first_diff)[["d_income"]], FE = b_within)       # T = 3: need not agree

# With T = 2 (periods 2 and 3 only), FD and FE coincide exactly
fe_t2 <- plm(liquor ~ income, data = filter(panel, year %in% 2:3),
             index = c("hh", "year"), model = "within")
fd_t2 <- lm(d_liquor ~ 0 + d_income, data = filter(panel_fd, year == 3))
c(FE_T2 = coef(fe_t2)[["income"]], FD_T2 = coef(fd_t2)[["d_income"]])


# Lab 1.4 · Are dummy variables and demeaning different estimators? -----------

fixed_dummies <- lm(liquor ~ income + factor(hh), data = panel)
coef(fixed_dummies)["income"]

fixed <- plm(liquor ~ income, data = panel, index = c("hh", "year"),
             model = "within", effect = "individual")
coef(fixed)

# [EXTRA] Frisch–Waugh–Lovell (POE5 §15.2.4, Exercise 15.16d):
# residualising on household dummies is the same as demeaning.
liquor_resid <- resid(lm(liquor ~ factor(hh), data = panel))
income_resid <- resid(lm(income ~ factor(hh), data = panel))
coef(lm(liquor_resid ~ 0 + income_resid))
all.equal(unname(income_resid), panel$income_within)

# [EXTRA] Same slope, different conventional standard errors
# (POE5 Examples 15.4–15.5, Exercise 15.16e). lm() on demeaned data uses
# NT - 1 degrees of freedom; demeaning has already used N household means,
# so the correct divisor is NT - N - K_S.
K_S <- 1
se_manual <- summary(within_manual)$coefficients["income_within", "Std. Error"]
se_lsdv   <- summary(fixed_dummies)$coefficients["income", "Std. Error"]
se_plm    <- summary(fixed)$coefficients["income", "Std. Error"]

c(manual           = se_manual,
  manual_corrected = se_manual * sqrt((NT - K_S) / (NT - N - K_S)),
  LSDV             = se_lsdv,
  plm              = se_plm)

# [EXTRA] Recovering household intercepts (reference material):
# alpha_i = mean(liquor_i) - b_FE * mean(income_i)
tibble(
  hh      = household_means$hh,
  by_hand = household_means$liquor_mean - b_within * household_means$income_mean,
  plm     = as.numeric(fixef(fixed))
) |>
  head()

# [EXTRA] F-test for household differences in intercepts (POE5 eq. 15.20,
# Exercise 15.16c). Assumes homoskedastic, uncorrelated errors.
anova(pooled, fixed_dummies)

# Interpretation: rejection says household intercepts differ. It does NOT say
# the pooled slope is biased. Bias requires the household effects to be
# correlated with income (Section 0). Heterogeneity is not endogeneity.


# Lab 1.5 · Does FE settle the substantive question? --------------------------
# Household-clustered CR1: G/(G-1) * (n-1)/(n-k), with k = 1 in the demeaned
# slope regression, and a t distribution with G - 1 = 39 degrees of freedom.

G <- N
vc_fe_cluster <- vcovCL(within_manual, cluster = ~hh, type = "HC1")

fe_inference <- coeftest(within_manual, vcov. = vc_fe_cluster, df = G - 1)
fe_inference
coefci(within_manual, vcov. = vc_fe_cluster, df = G - 1)

# The same numbers from the plm object ("sss" = Stata-style small-sample CR1)
vc_plm_cluster <- vcovHC(fixed, method = "arellano", type = "sss", cluster = "group")
coeftest(fixed, vcov. = vc_plm_cluster, df = G - 1)

# [EXTRA] Same estimator, different small-sample convention: clustering the
# LSDV regression counts the 40 dummies (k = 41) and inflates the SE.
# Do not mix conventions when comparing with other software.
sqrt(vcovCL(fixed_dummies, cluster = ~hh, type = "HC1")["income", "income"])

# [EXTRA] Alternative standard errors, pooled and FE (cf. POE5 Tables 15.3–15.4).
# Heteroskedasticity-robust (White) SEs are not valid for FE when T > 2
# (POE5 footnote 9); they are shown only for comparison.
tibble(
  model           = c("Pooled OLS", "Fixed effects"),
  estimate        = c(b_pooled, b_within),
  conventional    = c(sqrt(vcov(pooled)["income", "income"]), se_plm),
  heteroskedastic = c(sqrt(vcovHC(pooled, type = "HC1")["income", "income"]),
                      sqrt(vcovHC(within_manual, type = "HC1")[1, 1])),
  cluster         = c(sqrt(vcovCL(pooled, cluster = ~hh, type = "HC1")["income", "income"]),
                      sqrt(vc_fe_cluster[1, 1]))
)

# Clustering does not always increase the SE: compare the two rows.
# With G = 40 clusters, the number of households (not 120) is the effective
# sample size for inference.

# Deliberate misconception: "We included household fixed effects, so we have
# established a positive causal income effect." Name the two separate problems.


# 4. The AI audit -------------------------------------------------------------
# The numbers used in the audit prompt, reproduced from the objects above:
round(c(pooled          = b_pooled,
        between         = b_between,
        fixed_effects   = b_within,
        within_share    = within_share,
        clustered_se_fe = sqrt(vc_fe_cluster[1, 1])), 5)

# [EXTRA] Stress test in a simulation where the true effect is known.
# Scenario A: stable preferences are correlated with income.
# Scenario B: as A, plus a time-varying health shock that lowers income and
#             raises liquor expenditure.
# The true within-household income effect is beta = 0.02 in both scenarios.
set.seed(2026)
n_hh <- 20000  # try n_hh <- 40: noisier estimates; more households do not remove pooled bias
beta <- 0.02

sim <- expand_grid(hh = 1:n_hh, year = 1:3) |>
  group_by(hh) |>
  mutate(
    preference  = rnorm(1),                                # stable, unobserved
    income_base = 70 + 15 * preference + rnorm(1, sd = 10) # related to preference
  ) |>
  ungroup() |>
  mutate(
    health_shock = rnorm(n()),                             # time-varying, unobserved
    income_A = income_base + rnorm(n(), sd = 4),
    liquor_A = 1 + preference + beta * income_A + rnorm(n(), sd = 0.5),
    income_B = income_A - 3 * health_shock,
    liquor_B = 1 + preference + beta * income_B + 0.3 * health_shock +
               rnorm(n(), sd = 0.5)
  )

c(true      = beta,
  pooled_A  = coef(lm(liquor_A ~ income_A, data = sim))[["income_A"]],
  FE_A      = coef(plm(liquor_A ~ income_A, data = sim,
                       index = c("hh", "year"), model = "within"))[["income_A"]],
  pooled_B  = coef(lm(liquor_B ~ income_B, data = sim))[["income_B"]],
  FE_B      = coef(plm(liquor_B ~ income_B, data = sim,
                       index = c("hh", "year"), model = "within"))[["income_B"]])

# FE removes the stable confounder (A) but not the time-varying one (B).


# 5. Exit problem -------------------------------------------------------------
# A time-invariant household characteristic in the liquor data:
panel <- panel |>
  mutate(high_income_hh = as.numeric(income_mean > median(income_mean)))

# plm silently drops it: only the income coefficient is reported
coef(plm(liquor ~ income + high_income_hh, data = panel,
         index = c("hh", "year"), model = "within"))

# lm() with dummies drops whichever collinear column comes LAST.
# Here a household dummy is dropped, and high_income_hh gets a number
# that is not an effect of anything: it depends on the column order.
coef(lm(liquor ~ income + high_income_hh + factor(hh), data = panel))[1:3]
tail(coef(lm(liquor ~ income + factor(hh) + high_income_hh, data = panel)), 2)

# A dropped (or arbitrary) coefficient is not a zero effect: the variable
# has no within-household variation, so its separate effect is unidentified.


# Next: panel_data_02 ---------------------------------------------------------
# [EXTRA] Bridge (Mundlak, POE5 §15.4.3): adding the household mean of income
# to a pooled regression reproduces the FE slope exactly (balanced panel),
# and the coefficient on income_mean equals b_between - b_within.
mundlak <- lm(liquor ~ income + income_mean, data = panel)
coef(mundlak)
c(FE = b_within, between_minus_within = b_between - b_within)
