# PANEL DATA: AN ANNOTATED R CODING RESOURCE
# =============================================================================
# Consolidated from the Wooldridge scripts.
# Run this file section by section in RStudio (Ctrl+Enter).
# Sections depend on objects created earlier. No working-directory changes or
# workspace clearing are needed. Internet access is needed for the data URLs.
#
# ROUTE THROUGH THE FILE
#  1. Packages and panel structure
#  2. Pooled OLS: chemical firms
#  3. First differences with two periods
#  4. Fixed effects: within transformation and dummy variables
#  5. Clustered inference and random effects
#  6. Correlated random effects (Mundlak) and within-between models
#  7. Wage application: nonlinear terms, time effects and diagnostics
#  8. Transaction-price application: FE, RE and percentage interpretation
#  9. Crime application: controls and two-way fixed effects
# 10. Extensions: repeated cross sections and difference-in-differences
# 11. Extension: different firm slopes and a link to SUR
# 12. Advanced: pooled IV and fixed-effects IV for crime
# 13. Advanced: Hausman-Taylor estimation
#
# Sections 1-7 form the main teaching sequence. Sections 8-13 provide further
# applications and advanced material; this is a resource bank, not a requirement
# to cover everything in one sitting.
#
# MODEL TO KEEP IN MIND
# y_it = alpha + beta'x_it + u_i + e_it
# u_i: an unobserved unit characteristic, constant across that unit's observations.
# e_it: an observation-specific error, which can vary over time.
# Both FE and RE have a time-constant u_i. Conventional RE additionally assumes
# that u_i is uncorrelated with the entire regressor history. FE permits this
# correlation but still needs appropriate exogeneity of the time-varying error.
# Clustering changes inference; it cannot repair an inconsistent slope estimate.

# 1. Packages and panel structure ----
# Install once if needed; then leave this line commented out.
# install.packages(c("tidyverse", "broom", "plm", "lmtest",
#                    "sandwich", "car", "foreign", "systemfit"))
library(tidyverse)
library(broom)
library(plm)
library(lmtest)
library(sandwich)
library(car)

# Direct data loading, as in the original scripts.
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/nls_panel.rdata"))
head(nls_panel, 10)
glimpse(nls_panel)

# Each row is a person-year. The two indexes identify who and when.
nls_p <- pdata.frame(nls_panel, index = c("id", "year"))
pdim(nls_p)
is.pbalanced(nls_p)
# A balanced panel need not have equally spaced calendar observations.
sort(unique(nls_panel$year))

nls_panel %>%
  filter(id %in% head(sort(unique(nls_panel$id)), 3)) %>%
  select(id, year, lwage, educ, black, south, union, exper, tenure) %>%
  arrange(id, year)

# Check which regressors actually change within people. A variable can change
# for some people but remain constant for others (south and union, for example).
nls_panel %>%
  group_by(id) %>%
  summarise(across(c(educ, black, south, union, exper, tenure), n_distinct),
            .groups = "drop") %>%
  summarise(across(-id, ~ sum(.x > 1)))

# Use ordinary data frames for dplyr preparation, then declare the panel.
# On a pseries, use plm::lag() explicitly to avoid dplyr's row-based lag.
# Calendar gaps matter: plm::lag(x) means the previous indexed period, not
# automatically the previous available interview.

# 2. Pooled OLS: production in chemical firms ----
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/chemical2.rdata"))
chem_p <- pdata.frame(chemical2, index = c("firm", "year"))
pdim(chem_p)
head(chemical2)

# lsales, lcapital and llabor are natural logarithms.
# The capital coefficient is an elasticity, conditional on labor.
chem_ols <- lm(lsales ~ lcapital + llabor, data = chemical2)
chem_pool <- plm(lsales ~ lcapital + llabor, data = chem_p, model = "pooling")
summary(chem_ols)
all.equal(coef(chem_ols), coef(chem_pool))

# Pooling uses both differences across firms and changes within firms.
# If persistent productivity u_i is correlated with inputs, pooled slopes can
# mix input responses with omitted productivity differences.
100 * ((1.10)^coef(chem_ols)["lcapital"] - 1)
# Predicted percentage sales difference for 10% more capital, holding labor fixed.

# 3. First differences: start with two periods ----
chem_two <- chemical2 %>%
  filter(year %in% c(2005, 2006)) %>%
  arrange(firm, year)
chem_two_p <- pdata.frame(chem_two, index = c("firm", "year"))
pdim(chem_two_p)

# Delta y_it = beta' Delta x_it + Delta e_it: u_i - u_i = 0.
# diff() dispatches to the panel-aware method because these are pseries.
chem_two_p$d_sales <- diff(chem_two_p$lsales)
chem_two_p$d_capital <- diff(chem_two_p$lcapital)
chem_two_p$d_labor <- diff(chem_two_p$llabor)
head(chem_two_p[, c("firm", "year", "lsales", "d_sales", "d_capital")], 8)

chem_fd_manual <- lm(d_sales ~ 0 + d_capital + d_labor, data = chem_two_p)
chem_fd <- plm(lsales ~ 0 + lcapital + llabor,
               data = chem_two_p, model = "fd")
summary(chem_fd_manual)
summary(chem_fd)

# Why no intercept? We are differencing a levels model WITHOUT time effects.
# A common change between the two dates WOULD give an intercept in the FD
# regression. Adding one is a modeling choice, not an algebraic error.
chem_fd_time <- lm(d_sales ~ d_capital + d_labor, data = chem_two_p)
summary(chem_fd_time)

# 4. Fixed effects: three ways to obtain the same slopes ----
# Within transformation: y_it - mean_i(y) = beta'(x_it - mean_i(x))
#                                      + e_it - mean_i(e).
# Again, the unit effect disappears. With T=2 and matching time-effect choices,
# FD and FE give identical slopes.
chem_fe_two <- plm(lsales ~ lcapital + llabor,
                  data = chem_two_p, model = "within")
cbind(FD = coef(chem_fd), FE = coef(chem_fe_two))
all.equal(unname(coef(chem_fd)), unname(coef(chem_fe_two)))

# Now use all three years. FE and FD need not coincide when T > 2.
chem_fe <- plm(lsales ~ lcapital + llabor, data = chem_p, model = "within")
summary(chem_fe)

# Make demeaning visible using ordinary dplyr operations.
chem_demeaned <- chemical2 %>%
  group_by(firm) %>%
  mutate(sales_w = lsales - mean(lsales),
         capital_w = lcapital - mean(lcapital),
         labor_w = llabor - mean(llabor)) %>%
  ungroup()
chem_manual <- lm(sales_w ~ 0 + capital_w + labor_w, data = chem_demeaned)

# LSDV = least squares dummy variables. One intercept per firm, common slopes.
chem_lsdv <- lm(lsales ~ 0 + factor(firm) + lcapital + llabor, data = chemical2)
cbind(Within = coef(chem_fe),
      Demeaned_OLS = unname(coef(chem_manual)),
      Dummy_variables = coef(chem_lsdv)[c("lcapital", "llabor")])
# The manual demeaned lm gives the slopes, but its conventional degrees of
# freedom do not account for estimated firm effects. Use plm for FE inference.

# fixef(type="level") returns unit intercepts, not necessarily zero-mean effects.
head(fixef(chem_fe, type = "level"))
summary(fixef(chem_fe, type = "level"))
ggplot(data.frame(intercept = as.numeric(fixef(chem_fe, type = "level"))),
       aes(intercept)) +
  geom_histogram(bins = 25, fill = "steelblue", color = "white") +
  labs(title = "Estimated firm intercepts", x = "Firm intercept", y = "Count")

# Classical test: are all firm intercepts equal? This is not an FE-versus-RE test.
pFtest(chem_fe, chem_pool)
# Show the same F statistic without hard-coded sample sizes.
sse_pool <- sum(residuals(chem_pool)^2)
sse_fe <- sum(residuals(chem_fe)^2)
n_firms <- n_distinct(chemical2$firm)
f_manual <- ((sse_pool - sse_fe) / (n_firms - 1)) /
  (sse_fe / df.residual(chem_fe))
f_manual
pf(f_manual, n_firms - 1, df.residual(chem_fe), lower.tail = FALSE)
# This classical F reference distribution uses the usual spherical-error
# assumptions. Rejection indicates unequal intercepts; it does not prove
# regressor exogeneity or establish causality.

# 5. Clustered inference and random effects ----
# Retain chemical3 because the original examples use it for this comparison.
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/chemical3.rdata"))
chem3_p <- pdata.frame(chemical3, index = c("firm", "year"))
pdim(chem3_p)
chem3_pool <- plm(lsales ~ lcapital + llabor, data = chem3_p, model = "pooling")
chem3_fe <- plm(lsales ~ lcapital + llabor, data = chem3_p, model = "within")

# Arellano covariance: heteroskedasticity and arbitrary dependence WITHIN firms.
# The independent sampling units for this inference are firms.
# white1/white2 do not provide the same within-firm correlation allowance.
chem3_pool_v <- plm::vcovHC(chem3_pool, method = "arellano",
                           type = "HC1", cluster = "group")
chem3_fe_v <- plm::vcovHC(chem3_fe, method = "arellano",
                         type = "HC1", cluster = "group")
coeftest(chem3_pool)
coeftest(chem3_pool, vcov. = chem3_pool_v)
coeftest(chem3_fe, vcov. = chem3_fe_v)
# Notice: changing the covariance matrix changes SEs, not coefficients.
# For a plain lm, sandwich::vcovCL(model, cluster = data$firm) is the
# corresponding explicit clustering route; cluster="group" alone is not enough.

chem3_re <- plm(lsales ~ lcapital + llabor, data = chem3_p,
                model = "random", random.method = "swar")
summary(chem3_re)
ercomp(chem3_re)
# RE quasi-demeans: y_it - theta*ybar_i, and likewise for the regressors.
# theta=0 gives pooling; theta approaching 1 approaches the within transformation.
# RE retains (1-theta)*u_i, so correlation between u_i and regressors matters.
# Random effects does NOT mean a different u_i is drawn each year.

cbind(Pooled = coef(chem3_pool)[c("lcapital", "llabor")],
      FE = coef(chem3_fe), RE = coef(chem3_re)[c("lcapital", "llabor")])

# LM: null is zero individual-effect variance under the RE error structure.
# This is NOT a general heteroskedasticity test.
plmtest(chem3_pool, effect = "individual", type = "bp")
# Classical Hausman comparison: FE and RE slopes should agree under the null
# and maintained model assumptions. Non-rejection does not prove RE exogeneity.
phtest(chem3_fe, chem3_re)

# 6. Correlated random effects: the Mundlak approach ----
# Model u_i = gamma' xbar_i + a_i, with a_i orthogonal to included regressors.
# Add firm means to RE. This makes the permitted correlation explicit and
# gives a simple alternative to the panelr implementation in older scripts.
chem3_cre_data <- chemical3 %>%
  group_by(firm) %>%
  mutate(capital_mean = mean(lcapital), labor_mean = mean(llabor),
         capital_within = lcapital - capital_mean,
         labor_within = llabor - labor_mean) %>%
  ungroup()
chem3_cre_p <- pdata.frame(chem3_cre_data, index = c("firm", "year"))
chem3_cre <- plm(lsales ~ lcapital + llabor + capital_mean + labor_mean,
                 data = chem3_cre_p, model = "random")
chem3_cre_v <- plm::vcovHC(chem3_cre, method = "arellano",
                          type = "HC1", cluster = "group")
coeftest(chem3_cre, vcov. = chem3_cre_v)
cbind(FE = coef(chem3_fe), CRE = coef(chem3_cre)[c("lcapital", "llabor")])

# Joint test of mean terms: evidence against the conventional RE restriction.
linearHypothesis(chem3_cre, c("capital_mean = 0", "labor_mean = 0"),
                 vcov. = chem3_cre_v, test = "Chisq")

chem3_wb <- plm(lsales ~ capital_within + labor_within + capital_mean + labor_mean,
                data = chem3_cre_p, model = "random")
summary(chem3_wb)
# In this parameterization, within coefficients describe changes within firms;
# mean coefficients describe differences across firms. In the CRE version,
# beta_between = beta_within + gamma. These are NOT automatically short-run
# and long-run effects: no dynamic adjustment model has been estimated.

# 7. Wage application: nonlinear terms and time effects ----
# educ: schooling; exper: work experience; tenure: current-job tenure;
# south and union: indicator variables; black: the dataset's race indicator.
# Use a common complete sample for comparisons and construct squares first.
wage_data <- nls_panel %>%
  select(id, year, lwage, educ, exper, tenure, black, south, union) %>%
  drop_na() %>%
  mutate(exper2 = exper^2, tenure2 = tenure^2)
wage_p <- pdata.frame(wage_data, index = c("id", "year"))
wage_pool <- plm(lwage ~ educ + exper + exper2 + tenure + tenure2 + black + south + union,
                 data = wage_p, model = "pooling")
wage_fe <- plm(lwage ~ exper + exper2 + tenure + tenure2 + south + union,
               data = wage_p, model = "within")
wage_re <- plm(lwage ~ educ + exper + exper2 + tenure + tenure2 + black + south + union,
               data = wage_p, model = "random", random.method = "swar")

# Schooling and race do not change within these individuals; their separate
# slopes cannot be recovered by FE. This does NOT imply that their effects are zero.
wage_fe_v <- plm::vcovHC(wage_fe, method = "arellano", type = "HC1", cluster = "group")
coeftest(wage_fe, vcov. = wage_fe_v)
coeftest(wage_pool, vcov. = plm::vcovHC(wage_pool, method = "arellano",
                                      type = "HC1", cluster = "group"))
coeftest(wage_re, vcov. = plm::vcovHC(wage_re, method = "arellano",
                                    type = "HC1", cluster = "group"))
phtest(wage_fe, wage_re)

# Interpret coefficients by NAME, not their positions in an output table.
100 * (exp(coef(wage_fe)["union"]) - 1)
# Exact proportional contrast on the exponentiated log scale for union 0 -> 1.
# Without additional retransformation assumptions this is not automatically
# the same percentage contrast in the arithmetic conditional mean wage.

# With a quadratic, the slope depends on experience.
# d log(wage) / d exper = b_exper + 2*b_exper2*exper.
100 * (coef(wage_fe)["exper"] + 2 * coef(wage_fe)["exper2"] * 10)
# Approximate percentage slope at ten years of experience.

# Time effects absorb common shocks in each observed year.
wage_twfe <- plm(lwage ~ exper + exper2 + tenure + tenure2 + south + union,
                 data = wage_p, model = "within", effect = "twoways")
coeftest(wage_twfe, vcov. = plm::vcovHC(wage_twfe, method = "arellano",
                                      type = "HC1", cluster = "group"))
# An equivalent way is individual FE plus factor(year).
wage_year_dummies <- plm(lwage ~ exper + exper2 + tenure + tenure2 + south + union + factor(year),
                        data = wage_p, model = "within")
pFtest(wage_year_dummies, wage_fe)
# The F test is classical; shared year effects are substantively useful even
# when a particular test lacks power. Exact collinearity can drop regressors.

# Serial dependence and cross-person dependence are separate questions.
pwartest(wage_fe)
pcdtest(wage_year_dummies, test = "cd")
# With irregular interviews, interpret serial dependence over the observed waves
# carefully. Person clustering allows within-person dependence regardless of gaps.
# Person clustering does not allow arbitrary correlation ACROSS people.
# With only five waves, do not treat Driscoll-Kraay as an automatic solution.

# Within and between PLOTS: demean BOTH axes for the within relationship.
wage_plot <- wage_data %>%
  group_by(id) %>%
  mutate(lwage_w = lwage - mean(lwage), exper_w = exper - mean(exper)) %>%
  ungroup()
ggplot(wage_plot, aes(exper_w, lwage_w)) +
  geom_point(alpha = 0.12) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Within-person variation", x = "Experience minus person mean",
       y = "Log wage minus person mean")
wage_means <- wage_data %>%
  group_by(id) %>%
  summarise(exper = mean(exper), lwage = mean(lwage), .groups = "drop")
ggplot(wage_means, aes(exper, lwage)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Between-person variation", x = "Mean experience", y = "Mean log wage")
# These are descriptive bivariate plots, not the partial slopes from the full model.

# Wage first-difference exercise retained from the original scripts.
# The final two observed years are adjacent in this dataset.
wage_two <- wage_data %>% filter(year %in% tail(sort(unique(year)), 2))
wage_two_p <- pdata.frame(wage_two, index = c("id", "year"))
wage_fd <- plm(lwage ~ 0 + exper, data = wage_two_p, model = "fd")
summary(wage_fd)
# This intentionally simple model omits other changing wage determinants.

# 8. Transaction prices: a second FE/RE application ----
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mexican.rdata"))
# Repeated transactions for sex workers in Mexico. id identifies the worker;
# trans indexes transactions, not calendar years.
# lnprice: log transaction price; nocondom: unprotected transaction;
# bar/street: meeting location (othersite is the omitted category);
# rich/regular/alcohol: client characteristics;
# attractive/school/age: worker characteristics as recorded in the dataset.
mex_data <- mexican %>%
  select(id, trans, lnprice, bar, street, nocondom, rich, regular, alcohol,
         attractive, school, age) %>% drop_na()
mex_p <- pdata.frame(mex_data, index = c("id", "trans"))
pdim(mex_p)
mexican %>% summarise(across(c(bar, street, othersite), mean))

mex_pool <- plm(lnprice ~ bar + street + nocondom + rich + regular + alcohol + attractive + school + age,
                data = mex_p, model = "pooling")
mex_fe <- plm(lnprice ~ bar + street + nocondom + rich + regular + alcohol,
              data = mex_p, model = "within")
mex_re <- plm(lnprice ~ bar + street + nocondom + rich + regular + alcohol + attractive + school + age,
              data = mex_p, model = "random")
mex_fe_v <- plm::vcovHC(mex_fe, method = "arellano", type = "HC1", cluster = "group")
coeftest(mex_fe, vcov. = mex_fe_v)
phtest(mex_fe, mex_re)

# Match common slopes by name. Never subtract hard-coded column ranges.
mex_terms <- names(coef(mex_fe))
cbind(Pooled = coef(mex_pool)[mex_terms], FE = coef(mex_fe), RE = coef(mex_re)[mex_terms])
100 * (exp(coef(mex_fe)) - 1)

# Clustered, large-sample 95% interval, transformed to the percentage scale.
mex_b <- coef(mex_fe)["nocondom"]
mex_se <- sqrt(mex_fe_v["nocondom", "nocondom"])
mex_ci <- mex_b + qnorm(c(0.025, 0.975)) * mex_se
100 * (exp(c(estimate = mex_b, lower = mex_ci[1], upper = mex_ci[2])) - 1)
# Do not copy significance conclusions from older output: inspect this interval.
# FE controls stable worker differences, but transaction-specific bargaining,
# risk and selection can still confound a causal interpretation.

# 9. Crime: additional controls and two-way FE ----
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/crime.rdata"))
crime_p <- pdata.frame(crime, index = c("county", "year"))
pdim(crime_p)
# Log variables: crime rate, arrest/conviction/imprisonment probabilities,
# average sentence length, and manufacturing wages. Slopes are elasticities.
crime_pool <- plm(lcrmrte ~ lprbarr + lprbconv + lprbpris + lavgsen + lwmfg,
                  data = crime_p, model = "pooling")
crime_fe <- plm(lcrmrte ~ lprbarr + lprbconv + lprbpris + lavgsen + lwmfg,
                data = crime_p, model = "within")
coeftest(crime_fe, vcov. = plm::vcovHC(crime_fe, method = "arellano",
                                     type = "HC1", cluster = "group"))
pFtest(crime_fe, crime_pool)
crime_fe_controls <- update(crime_fe, . ~ . + ldensity + lpctymle)
crime_twfe <- update(crime_fe_controls, effect = "twoways")
coeftest(crime_twfe, vcov. = plm::vcovHC(crime_twfe, method = "arellano",
                                       type = "HC1", cluster = "group"))
# FE does not solve reverse causality: crime can itself influence policing.
# A deterrence elasticity alone is also insufficient for a cost-benefit ranking.

# Show lag, difference, between and within operations using the Wooldridge data.
# Read the Stata-data.
library(foreign)
crime4 <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/crime4.dta")
crime4_p <- pdata.frame(crime4, index = c("county", "year"))
crime4_p$crime_lag <- plm::lag(crime4_p$crmrte)
crime4_p$crime_diff <- diff(crime4_p$crmrte)
crime4_p$crime_between <- plm::Between(crime4_p$crmrte)
crime4_p$crime_within <- plm::Within(crime4_p$crmrte)
head(crime4_p[, c("county", "year", "crmrte", "crime_lag", "crime_diff",
                  "crime_between", "crime_within")], 14)

# Two-period unemployment/crime exercise. The original data are ordered as
# consecutive two-row blocks for 46 cities; pdata.frame(index=46) uses that
# documented ordering to create unit and time indexes. Do not reorder first.
crime2 <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/crime2.dta")
head(crime2)
table(crime2$year)
crime2_p <- pdata.frame(crime2, index = 46)
crime2_p$d_crime <- diff(crime2_p$crmrte)
crime2_p$d_unem <- diff(crime2_p$unem)
crime2_fd_manual <- lm(d_crime ~ d_unem, data = crime2_p)
crime2_fd <- plm(crmrte ~ unem, data = crime2_p, model = "fd")
cbind(Manual = coef(crime2_fd_manual), plm = coef(crime2_fd))
# Here the FD intercept is retained: it allows a common change in crime.

# 10. Repeated cross sections and difference-in-differences ----
# These examples are related extensions: they do not track the same people or
# houses over time, so do not apply person/house FE to invented identifiers.
cps <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/cps78_85.dta")
cps_model <- lm(lwage ~ y85 * (educ + female) + exper + I(exper^2 / 100) + union,
                data = cps)
coeftest(cps_model, vcov. = sandwich::vcovHC(cps_model, type = "HC1"))
# y85:educ measures the change in the schooling slope between survey years.
# y85:female measures the change in the conditional female wage differential.

kielmc <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/kielmc.dta")
# Housing prices near a proposed incinerator, before and after the event.
# nearinc: nearby location; y81: post-period indicator; rprice: real house price.
kielmc %>% group_by(year, nearinc) %>%
  summarise(mean_price = mean(rprice), n = n(), .groups = "drop")
did_basic <- lm(rprice ~ nearinc * y81, data = kielmc)
did_age <- update(did_basic, . ~ . + age + I(age^2))
did_controls <- update(did_age, . ~ . + intst + land + area + rooms + baths)
summary(did_basic)
summary(did_controls)
coef(did_basic)["nearinc:y81"]
# Interaction = (near after - near before) - (far after - far before).
# Causal interpretation needs parallel counterfactual trends, no spillovers,
# and suitable control of changes in sample composition. One pre-period cannot
# establish parallel trends empirically. House-level HC errors do not capture
# all shared geographic shocks; these regressions illustrate the estimator.
did_log <- lm(log(rprice) ~ nearinc * y81, data = kielmc)
100 * (exp(coef(did_log)["nearinc:y81"]) - 1)
# Ratio-of-ratios on the exponentiated log-price scale.

# 11. Different slopes across firms, and a link to SUR ----
# Retains the distinct Grunfeld material from the older Chapter 15 scripts.
# Use the built-in plm data directly; the first two firms form the comparison.
data("Grunfeld", package = "plm")
grun_two <- Grunfeld %>% filter(firm %in% c(1, 2)) %>% mutate(firm = factor(firm))
grun_common <- lm(inv ~ value + capital, data = grun_two)
grun_intercepts <- lm(inv ~ firm + value + capital, data = grun_two)
grun_slopes <- lm(inv ~ firm * (value + capital), data = grun_two)
anova(grun_common, grun_intercepts, grun_slopes)
# Common slopes with different intercepts is FE. Interacting firm with inputs
# ALSO allows different slopes. These are distinct restrictions.
# This ANOVA is a classical comparison, not a serial-correlation-robust test.

# SUR allows contemporaneously correlated errors across the two firm equations.
# It does not by itself allow arbitrary serial correlation within each equation.
library(systemfit)
grun_wide <- grun_two %>%
  select(year, firm, inv, value, capital) %>%
  pivot_wider(names_from = firm, values_from = c(inv, value, capital))
grun_system <- list(firm1 = inv_1 ~ value_1 + capital_1,
                    firm2 = inv_2 ~ value_2 + capital_2)
grun_sur <- systemfit(grun_system, method = "SUR", data = as.data.frame(grun_wide))
summary(grun_sur)
linearHypothesis(grun_sur, c("firm1_value_1 = firm2_value_2",
                             "firm1_capital_1 = firm2_capital_2"))
# The equation labels above deliberately make cross-equation restrictions readable.

# 12. Advanced: fixed effects plus instrumental variables ----
# FE removes stable county characteristics. IV addresses correlation of selected
# regressors with the remaining error. These solve different problems.
# Candidate excluded instruments in the original exercise: ltaxpc and lmix.
# Their relevance and exclusion restrictions must be argued, not assumed true
# because the software accepts the formula.
crime_iv_data <- crime %>%
  select(county, year, lcrmrte, lpolpc, lprbarr, lprbconv, lavgsen, lwmfg,
         west, urban, ltaxpc, lmix) %>% drop_na() %>% arrange(county, year)
crime_iv_p <- pdata.frame(crime_iv_data, index = c("county", "year"))

# First treat police per capita (lpolpc) as the single endogenous regressor.
# Terms after | are instruments: include all exogenous regressors there too.
crime_iv_pool <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban |
                     lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                    data = crime_iv_p, model = "pooling")
crime_iv_fe <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg |
                   lprbarr + lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                  data = crime_iv_p, model = "within")
coeftest(crime_iv_pool, vcov. = plm::vcovHC(crime_iv_pool, method = "arellano",
                                          type = "HC1", cluster = "group"))
coeftest(crime_iv_fe, vcov. = plm::vcovHC(crime_iv_fe, method = "arellano",
                                        type = "HC1", cluster = "group"))
# Time-invariant west and urban disappear under within transformation.

# The relevant first stage must also remove county effects.
police_first <- plm(lpolpc ~ lprbarr + lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                    data = crime_iv_p, model = "within")
linearHypothesis(police_first, c("ltaxpc = 0", "lmix = 0"), test = "Chisq",
                 vcov. = plm::vcovHC(police_first, method = "arellano",
                                     type = "HC1", cluster = "group"))
# Joint significance is evidence of relevance, not proof of strong identification
# or of instrument validity. Do not compare this chi-square with an F>10 rule.

# Residual-inclusion exogeneity check on the SAME ordered complete sample.
crime_iv_p$police_residual <- as.numeric(residuals(police_first))
crime_control_function <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + police_residual,
                              data = crime_iv_p, model = "within")
linearHypothesis(crime_control_function, "police_residual = 0", test = "Chisq",
                 vcov. = plm::vcovHC(crime_control_function, method = "arellano",
                                     type = "HC1", cluster = "group"))
# This is an augmented FE regression, not an IV regression that instruments
# the first-stage residual. Interpret it under valid-instrument assumptions.

# Now also allow arrest probability to be endogenous.
crime_iv_two <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg |
                    lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                   data = crime_iv_p, model = "within")
police_first_two <- plm(lpolpc ~ lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                        data = crime_iv_p, model = "within")
arrest_first_two <- plm(lprbarr ~ lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                        data = crime_iv_p, model = "within")
summary(police_first_two)
summary(arrest_first_two)
coeftest(crime_iv_two, vcov. = plm::vcovHC(crime_iv_two, method = "arellano",
                                         type = "HC1", cluster = "group"))
# Two excluded instruments for two endogenous regressors: exactly identified
# if the rank condition holds. No overidentification test is available here.
# Separate significant first stages do not establish joint identification strength.
# The old hand-calculated Sargan N*R^2 test is omitted: its classical reference
# distribution does not generally survive within-county error dependence.

# 13. Advanced: Hausman-Taylor (HT) ----
# HT uses internal instruments constructed from the panel regressors.
# It can identify coefficients on some time-invariant variables correlated with
# u_i, if enough valid time-varying exogenous variables provide identification.
# ALL regressors must still be suitably exogenous to the idiosyncratic error.
# This is a stronger structured argument than simply choosing FE or RE.
#
# Current plm syntax: three RHS parts.
# 1: all regressors; 2: regressors uncorrelated with u_i;
# 3: time-varying regressors allowed to correlate with u_i.
# Source: https://search.r-project.org/CRAN/refmans/plm/html/pht.html

# Wage classification retained from the older exercise:
# educ (time-invariant) and south (time-varying) may correlate with u_i;
# exper, exper2, tenure, tenure2, union and black are assumed uncorrelated with u_i.
wage_ht <- plm(lwage ~ educ + exper + exper2 + tenure + tenure2 + black + south + union |
               exper + exper2 + tenure + tenure2 + union + black | south,
              data = wage_p, model = "random", random.method = "ht",
              inst.method = "baltagi")
summary(wage_ht)
# The credibility of that classification is the substantive question.

# Transaction-price classification from the original exercise:
# nocondom may correlate with stable worker traits, but is assumed exogenous
# to transaction-specific errors. All other included regressors are placed
# in the exogenous group; this restriction is open to criticism.
mex_ht <- plm(lnprice ~ bar + street + nocondom + rich + regular + alcohol + attractive + school + age |
              bar + street + rich + regular + alcohol + attractive + school + age | nocondom,
             data = mex_p, model = "random", random.method = "ht",
             inst.method = "baltagi")
summary(mex_ht)
100 * (exp(sum(coef(mex_ht)[c("nocondom", "attractive", "school")])) - 1)
# Combined contrast requires summing log coefficients before exponentiating.
# It does not establish a causal price premium without the model assumptions.

# QUESTIONS TO ANSWER AFTER RUNNING THE CODE ----
# 1. Which variation identifies the pooled, FE and between coefficients?
# 2. Why do FD and FE coincide for two periods under matching specifications?
# 3. What does clustering change, and what does it leave unchanged?
# 4. What extra restriction does conventional RE impose relative to FE?
# 5. What do the Mundlak mean terms represent? What does their joint test assess?
# 6. Why do education and race disappear in the individual FE wage model?
# 7. Can FE remove reverse causality between crime and police resources?
# 8. What must be true of the instruments before the IV slopes are credible?
# 9. Which HT variable classifications would you find difficult to defend?
#
# Consolidation notes:
# Repeated estimation/printing blocks and repeated data downloads were merged.
# Obsolete plm.data(), do(), machine-specific setwd(), workspace clearing,
# silent try() wrappers, and orphaned mlogit/Train examples were removed.
# Within-between material is expressed directly using plm and explicit means.
# Empirical conclusions are left to the displayed estimates and suitable SEs.
