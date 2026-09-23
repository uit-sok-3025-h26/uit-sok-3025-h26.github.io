#' ---
#' title: "Panel Data Econometrics in R: from pooled cross sections to panel IV"
#' subtitle: "R coding session following the two panel data lectures"
#' ---
#'
#' This script is a standalone resource. Run it from top to bottom.
#' The main reference is Principles of Econometrics, 5th ed. (POE5), Chapter 15.
#' Some examples come from POE4 Chapter 15, from Wooldridge's
#' Introductory Econometrics (Chapter 13), and from Heiss' "Using R for
#' Introductory Econometrics" (URFIE).
#'
#' The script moves from simple to more complicated models:
#'
#'  Part 1.  Pooled cross sections over time (OLS with time dummies and interactions)
#'  Part 2.  The structure of panel data, and panel data operations (lag, diff, between, within)
#'  Part 3.  The first-difference (FD) estimator
#'  Part 4.  The fixed effects (FE) / within estimator, and LSDV
#'  Part 5.  Testing for fixed effects (pooled OLS vs. FE)
#'  Part 6.  Cluster-robust standard errors
#'  Part 7.  The random effects (RE) estimator
#'  Part 8.  The Hausman test (FE vs. RE)
#'  Part 9.  A complete workflow: the wage equation, with diagnostics
#'  Part 10. Within and between effects (within-between model and plots)
#'  Part 11. The Hausman-Taylor estimator
#'  Part 12. Sets of regression equations (SUR)
#'  Part 13. Applied exercises: POE5 15.18, 15.22, 15.29 and 15.30 (panel IV)

rm(list = ls())

#' Load packages. Note that plm is loaded LAST on purpose:
#' plm has its own lag(), lead() and between() functions for panel data,
#' and they are masked by dplyr if dplyr is loaded after plm.
#' To be safe, we also write plm::lag() explicitly below.
suppressPackageStartupMessages(require(pacman))
p_load(knitr, broom, dplyr, mosaic, AER, lmtest, car, stargazer, foreign,
       ggplot2, panelr, sem, systemfit, plm)


# =============================================================================
#' # Part 1. Pooled cross sections over time
# =============================================================================
#' Before we look at "true" panel data, where we follow the SAME individuals
#' over time, we look at independent cross sections sampled at different
#' points in time. Different individuals in each period, so no individual
#' effects. We can still use time dummies and interactions to study how
#' relationships change over time.

#' Wooldridge data sets are available here:
# browseURL("http://fmwww.bc.edu/ec-p/data/wooldridge/datasets.list.html")

# -----------------------------------------------------------------------------
#' Wooldridge, Example 13.2: Changes in the return to education and the gender wage gap.
#' CPS data for 1978 and 1985. y85 is a dummy for 1985.
#' The interaction y85*(educ+female) lets the return to education and the
#' gender gap differ between the two years.
cps <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/cps78_85.dta")

#' Detailed OLS results including interaction terms
summary(lm(lwage ~ y85*(educ+female) + exper + I((exper^2)/100) + union, data = cps))

# -----------------------------------------------------------------------------
#' Wooldridge, Example 13.3: Effect of a garbage incinerator's location on house prices.
#' This is a difference-in-differences (DiD) design with two cross sections.
kielmc <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/kielmc.dta")
head(kielmc)

#' Separate regressions for 1978 and 1981
table(kielmc$year)
kielmc %>% group_by(year) %>% do(tidy(lm(rprice ~ nearinc, data = .)))
kielmc %>% group_by(year) %>% do(glance(lm(rprice ~ nearinc, data = .))) %>% round(., 2)

#' Joint regression including an interaction term.
#' The coefficient on nearinc:y81 is the difference-in-differences estimate.
#' Table 13.2, model (1)
kielmc %>% do(tidy(lm(rprice ~ nearinc*y81, data = .)))
kielmc %>% do(glance(lm(rprice ~ nearinc*y81, data = .))) %>% round(., 2)

#' Table 13.2, model (2)
kielmc %>% do(tidy(lm(rprice ~ nearinc*y81 + age + I(age^2), data = .)))
kielmc %>% do(glance(lm(rprice ~ nearinc*y81 + age + I(age^2), data = .))) %>% round(., 2)

#' Table 13.2, model (3)
names(kielmc)
kielmc %>% do(tidy(lm(rprice ~ nearinc*y81 + age + I(age^2) + intst + land + area + rooms + baths, data = .)))
kielmc %>% do(glance(lm(rprice ~ nearinc*y81 + age + I(age^2) + intst + land + area + rooms + baths, data = .))) %>% round(., 2)

#' Log-lin model. The DiD coefficient is now (approximately) a percentage effect.
#' The exact percentage effect of a dummy in a log-lin model is 100*(exp(b)-1).
did.log <- lm(log(rprice) ~ nearinc*y81, data = kielmc)
tidy(did.log)
glance(did.log) %>% round(., 2)
100*(exp(coef(did.log)["nearinc:y81"]) - 1)


# =============================================================================
#' # Part 2. The structure of panel data
# =============================================================================

# -----------------------------------------------------------------------------
#' ## Example 15.1 (POE5)
#' This example lists observations on several variables in a microeconometric panel of individuals.
#' The nls_panel.rdata dataset includes a subset of National Longitudinal Survey which
#' is conducted by the U.S. Department of Labor. The database includes observations on women, who
#' in 1968, were between the ages of 14 and 24. It then follows them through time, recording various
#' aspects of their lives annually until 1973 and bi-annually afterwards. Here, we use a sub-sample of
#' N=716 women who were interviewed in 1982, 1983, 1985, 1987 and 1988. The sample consists of women
#' who were employed, and whose schooling was completed, when interviewed. The panel is balanced and
#' there are 3580 total observations.

#' Data definition file:
# browseURL("http://www.principlesofeconometrics.com/poe5/data/def/nls_panel.def")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/nls_panel.rdata"))
head(nls_panel, 10)
# View(nls_panel)

str(nls_panel)
glimpse(nls_panel)

#' Create a panel data frame using the plm panel data package.
#' To use panel data functions, we must first declare the data's structure.
#' 'id' is the individual (cross-sectional) identifier.
#' 'year' is the time-series identifier.
nlsPDF <- pdata.frame(nls_panel, index = c("id", "year"))

#' pdim() checks the panel dimensions and balance.
#' A 'balanced panel' means every individual is observed in every time period.
pdim(nlsPDF)

#' The structure now shows that variables are of class 'pseries',
#' which is a special panel-aware data type used by the 'plm' package.
glimpse(nlsPDF) # note pseries

#' Table 15.1: the first three individuals.
#' Notice how some variables (like educ, black, south) are time-invariant for these women,
#' while others (like lwage, exper, tenure) vary over time.
nlsPDF %>% dplyr::filter(as.integer(id) <= 3) %>%
  select(id, year, lwage, educ, south, black, union, exper, tenure)

#' Same table, formatted with kable()
nlsPDF %>% dplyr::filter(as.integer(id) <= 3) %>%
  select(id, year, lwage, educ, south, black, union, exper, tenure) %>%
  kable(., digits = 4, align = "c", caption = "Table 15.1: Panel data for three women")

# -----------------------------------------------------------------------------
#' ## Panel data calculations (URFIE, p. 202-203)
#' Once the data are a pdata.frame, R "knows" which observations belong to the same
#' individual. Lags, differences, individual means (between) and deviations from
#' individual means (within, "demeaning") are then computed within each individual.
#'
#' Note: plm::lag, plm::lead and plm::between are masked if dplyr and/or mosaic
#' are loaded after plm. We therefore write plm::lag() explicitly.
crime4 <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/crime4.dta")

#' Generate pdata.frame:
crime4.p <- pdata.frame(crime4, index = c("county", "year"))
pdim(crime4.p)

#' Calculations within the pdata.frame:
crime4.p$cr.l <- plm::lag(crime4.p$crmrte)   # lag, NA in the first year of each county
crime4.p$cr.d <- diff(crime4.p$crmrte)       # first difference
crime4.p$cr.B <- Between(crime4.p$crmrte)    # individual (county) mean, repeated each year
crime4.p$cr.W <- Within(crime4.p$crmrte)     # deviation from the county mean (demeaning)

#' Display selected variables for observations 1-16:
crime4.p[1:16, c("county", "year", "crmrte", "cr.l", "cr.d", "cr.B", "cr.W")]
#' Check: crmrte = cr.B + cr.W for every observation.

#' Note: diff() and plm::lag() do not work inside dplyr pipes on pseries (panel) objects:
# crime4.p %>% mutate(cr.d = diff(crmrte))   # does not do what you expect
#' Alternative: the panelr package supports mutate() with lags on panel data (see Part 3).


# =============================================================================
#' # Part 3. The first-difference (FD) estimator
# =============================================================================
#' The first set of models is based on the presence of fixed effects in the regression function as
#' shown in equation (15.9). When there are only two time periods, the data can be time-differenced
#' and OLS used to estimate the slopes of all time-varying regressors. Time-invariant variables and
#' the intercept drop out of the model upon differencing.
#'
#' Pedagogical notes on the difference estimator (2 periods):
#' - With logs, coefficients are elasticities of *changes*.
#' - Differencing removes time-invariant heterogeneity but can amplify measurement error.
#' - With exactly T=2, first differences and FE-within yield identical slope estimates (Part 4).

# -----------------------------------------------------------------------------
#' ## Two-period panel, Wooldridge crime2 data (First Difference Models)
crime2 <- read.dta("http://fmwww.bc.edu/ec-p/data/wooldridge/crime2.dta")
head(crime2)
names(crime2)

#' The data are stacked: 46 cities, two years. index=46 tells plm that there are
#' 46 individuals; plm then creates the variables 'id' and 'time'.
PDcrime2 <- pdata.frame(crime2, index = 46)
pdim(PDcrime2)

#' Manually calculate first differences:
PDcrime2$Dcrmrte <- diff(PDcrime2$crmrte)
PDcrime2$Dunem   <- diff(PDcrime2$unem)

#' Display selected variables for observations 1-6:
PDcrime2[1:6, c("id", "time", "year", "crmrte", "Dcrmrte", "unem", "Dunem")]

#' Cross-section regression in each year: unemployment seems to have no effect on crime
PDcrime2 %>% group_by(time) %>% do(tidy(lm(crmrte ~ unem, data = .)))

#' Estimate FD model with lm on differenced data:
tidy(lm(Dcrmrte ~ Dunem, data = PDcrime2))
glance(lm(Dcrmrte ~ Dunem, data = PDcrime2))

#' Estimate FD model with plm on the original data (same result):
tidy(plm(crmrte ~ unem, model = "fd", data = PDcrime2))

# -----------------------------------------------------------------------------
#' ## Example 15.2 The Difference Estimator (POE5)
#' This is illustrated in Example 15.2 in POE5. The data are included in the chemical2.rdata dataset.
#' This dataset contains sales, capital, and labor inputs for Chinese chemical firms. There are 3 years
#' of observations on 200 firms. The model to be estimated is a log-log model of sales.
# browseURL("http://www.principlesofeconometrics.com/poe5/data/def/chemical2.def")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/chemical2.rdata"))
head(chemical2, 10)

chemPDF <- pdata.frame(chemical2, index = c("firm", "year"))
pdim(chemPDF)

#' Pooled OLS on the two years 2005 and 2006 (ignores the firm effects):
chemical2 %>% filter(year %in% c(2005:2006)) %>% do(tidy(lm(lsales ~ lcapital + llabor, data = .)))

#' Create variables, use only 2005 and 2006 data. We call this panel chemPDF_T2 (T=2).
chemPDF_T2 <- pdata.frame(filter(chemical2, year %in% c(2005:2006)), index = c("firm", "year"))
pdim(chemPDF_T2)

#' Taking the time difference yields:
#' dlnsales = ln(sales_t) - ln(sales_{t-1}), which is (approximately) the growth rate of sales.
chemPDF_T2$dlnsales   <- diff(chemPDF_T2$lsales)
chemPDF_T2$dlncapital <- diff(chemPDF_T2$lcapital)
chemPDF_T2$dlnlabor   <- diff(chemPDF_T2$llabor)

#' The panelr package can use mutate() on lags (here dplyr::lag works within each firm):
chemical2 %>%
  panel_data(., id = firm, wave = year) %>%
  mutate(dlnsales   = lsales - dplyr::lag(lsales),
         dlncapital = lcapital - dplyr::lag(lcapital),
         dlnlabor   = llabor - dplyr::lag(llabor)) %>%
  select(year, firm, lsales, dlnsales, lcapital, dlncapital, llabor, dlnlabor) %>% head(., 10)

#' Equation in first differences of logs, no intercept.
#' Note the '0 +' syntax, which means we estimate without an intercept.
#' The intercept (and the fixed effect) is removed by differencing.
fd_model <- lm(dlnsales ~ 0 + dlncapital + dlnlabor, data = chemPDF_T2)
tidy(fd_model)
glance(fd_model)

#' Easier: let plm do the differencing for you
tidy(plm(lsales ~ 0 + lcapital + llabor, data = chemPDF_T2, model = "fd"))

# -----------------------------------------------------------------------------
#' ## Example 15.3 The Difference Estimator, wage equation (POE5)
#' The difference estimator is used to estimate a simple wage regression based on the nls_panel data.
#' Note that the nls_panel2 data does not exist in POE5.
#' We have to make a panel, using the last two years:
nls_panel %>% filter(year %in% c(87, 88)) -> nls_panel2

#' Create a panel data frame using plm panel data package
nls2PDF <- pdata.frame(nls_panel2, index = c("id", "year"))
pdim(nls2PDF)

#' Create variables for the first-differenced wage model
nls2PDF$dlnwage <- diff(nls2PDF$lwage)
nls2PDF$dexper  <- nls2PDF$exper - plm::lag(nls2PDF$exper)

#' Interpretation: the coefficient on dexper is a semi-elasticity of the change in ln(wage)
#' with respect to the change in experience (~ % change in wage per extra year).
summary(lm(dlnwage ~ 0 + dexper, data = nls2PDF))


# =============================================================================
#' # Part 4. The fixed effects (FE) / within estimator
# =============================================================================

# -----------------------------------------------------------------------------
#' ## Example 15.4, The within estimator (POE5)
#' In this example the within transformation is used to estimate a two-period fixed effects model
#' of Chinese chemical firms.
#' Remember (in log first diff):
summary(fd_model)

#' Note model now in logs, not in log diff. Notice model = "within".
example15_4 <- plm(lsales ~ lcapital + llabor, data = chemPDF_T2, model = "within")
summary(example15_4)
#' Same parameters as the difference estimator above, and the same standard errors.
#' The FE estimator de-means the data instead of differencing. With T=2 the two
#' transformations give identical slope estimates.

#' The 200 intercepts per firm (id)
summary(fixef(example15_4))

# -----------------------------------------------------------------------------
#' ## Example 15.5 T=3 (POE5)
#' In this example, the within transformation is used on the sample with T = 3 years of data.
#' With T>2, the FD and FE estimators are no longer identical, though both are consistent.
#' The FE estimator is more efficient if the errors e_it are not serially correlated.
example15_5 <- plm(lsales ~ lcapital + llabor, data = chemPDF, model = "within")
summary(example15_5)
summary(fixef(example15_5))

#' Interpretation (log-log model): a 10% increase in capital implies approximately a
#' 10*b percent increase in sales, where b is:
coef(example15_5)["lcapital"]

#' ## FE vs. LSDV
#' An equivalent way to estimate this model is using the least squares dummy variable
#' estimator (LSDV). Here, an indicator variable is created for each individual in
#' the sample. These are added to the model (dropping the intercept) and estimated by least squares.
#' There is a good reason why this formulation of the fixed effects model is not used more often:
#' it produces a ton of output. Since there are 200 firms in the data, 200 lines of extra output
#' would be sent to the screen. We therefore only print the slopes.
lsdv_chem <- lm(lsales ~ 0 + factor(firm) + lcapital + llabor, data = chemical2)
coef(lsdv_chem)[c("lcapital", "llabor")]

#' Quick check: equality of slopes from 'within' and LSDV
coef(example15_5) - coef(lsdv_chem)[c("lcapital", "llabor")]

#' Notes on FE vs. LSDV:
#' - summary(example15_5) reports the within R^2 (fit after demeaning by unit means).
#' - fixef(example15_5) gives the alpha_i (unit intercepts).
#' - LSDV (dummy per firm, no constant) yields identical slopes and fitted values.
#'
#' Notes on R^2 in FE models:
#' - Within: variation around each unit mean (most relevant for FE interpretation)
#' - Between: across-unit variation in unit means
#' - Overall: pooled variation ignoring the FE structure


# =============================================================================
#' # Part 5. Testing for fixed effects (pooled OLS vs. FE)
# =============================================================================

# -----------------------------------------------------------------------------
#' ## Example 15.6 (POE5)
#' Testing for unobserved heterogeneity, i.e., what is better, one intercept (for all)
#' or an intercept per id?
#'
#' H0: All firm-specific intercepts are equal (alpha_1 = alpha_2 = ... = alpha_N).
#'     (Pooled OLS is sufficient.)
#' H1: At least one intercept is different. (Fixed effects are necessary.)
pdim(chemPDF)

#' OLS pooled, one intercept, the restricted model
mod15_6OLS <- plm(lsales ~ lcapital + llabor, data = chemPDF, model = "pooling")
summary(mod15_6OLS)

#' Fixed effects, one intercept per id, the unrestricted model
mod15_6FE <- plm(lsales ~ lcapital + llabor, data = chemPDF, model = "within")
mean(fixef(mod15_6FE))

#' Do the F-test manually (equation 15.20)
glance(mod15_6OLS)
glance(mod15_6FE)

SSER <- sum(resid(mod15_6OLS)^2)   # restricted sum of squared residuals
SSEU <- sum(resid(mod15_6FE)^2)    # unrestricted sum of squared residuals
N  <- length(fixef(mod15_6FE))     # number of firms
NT <- nrow(chemPDF)                # number of observations
K  <- length(coef(mod15_6FE))      # number of slope parameters
#' F-value
((SSER - SSEU)/(N - 1)) / (SSEU/(NT - N - K))
#' Critical F
qf(0.99, N - 1, NT - N - K)
#' The hypothesis that the fixed effects are equal to one another is rejected at 1%.

#' The built-in test is simpler:
pFtest(mod15_6FE, mod15_6OLS)
#' Rejecting H0 means the alpha_i are not all equal, so FE is preferred over pooled OLS.


# =============================================================================
#' # Part 6. Cluster-robust standard errors
# =============================================================================

# -----------------------------------------------------------------------------
#' ## Example 15.7 & 15.8 Robust standard errors (POE5)
#' Robust covariances in panel data take into account the special nature of these data.
#' Specifically they account for autocorrelation within the observations on each individual and
#' they allow the variances for different individuals to vary.
#' Since panel data have both a time series and a cross-sectional dimension one might expect that,
#' in general, robust estimation of the covariance matrix would require handling both
#' heteroskedasticity and autocorrelation (at the same time).
# browseURL("http://www.principlesofeconometrics.com/poe5/data/def/chemical3.def")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/chemical3.rdata"))

chem3PDF <- pdata.frame(chemical3, index = c("firm", "year"))
pdim(chem3PDF)

#' Pooled OLS, estimated with plm so that R knows the panel structure
mod15_7OLS <- plm(lsales ~ lcapital + llabor, data = chem3PDF, model = "pooling")
summary(mod15_7OLS)

#' Correction using the built-in function in plm:
# vcovHC(x, method = c("arellano", "white1", "white2"),
#        type = c("HC0", "sss", "HC1", "HC2", "HC3", "HC4"),
#        cluster = c("group", "time"))
#' method = "white1"  : heteroskedasticity only (no within-firm correlation)
#' method = "arellano": heteroskedasticity + arbitrary serial correlation within firm (cluster-robust)
#' type = "sss"       : Stata-like small sample correction of the cluster-robust estimator
#'
#' WARNING: if you estimate the model with lm() instead of plm(), vcovHC() is the version
#' from the sandwich package. It does not know the panel structure and silently
#' ignores cluster = "group". Always use plm() when you want panel-robust standard errors.

sqrt(diag(vcov(mod15_7OLS)))                                                            # conventional
sqrt(diag(vcovHC(mod15_7OLS, method = "white1",   type = "HC0")))                       # heteroskedasticity only
sqrt(diag(vcovHC(mod15_7OLS, method = "arellano", type = "HC0", cluster = "group")))    # cluster-robust
sqrt(diag(vcovHC(mod15_7OLS, method = "arellano", type = "sss", cluster = "group")))    # cluster-robust, small sample

#' The full table with cluster-robust standard errors, t-values and p-values:
coeftest(mod15_7OLS, vcov = vcovHC(mod15_7OLS, method = "arellano", type = "HC0", cluster = "group"))

#' Using the fixed effects estimator
mod15_8FE <- plm(lsales ~ lcapital + llabor, data = chem3PDF, model = "within")
summary(mod15_8FE)
coeftest(mod15_8FE, vcov = vcovHC(mod15_8FE, method = "arellano", type = "HC0", cluster = "group"))

#' If T is moderate (e.g., >= 5) and there may be cross-sectional dependence,
#' consider Driscoll-Kraay standard errors:
# coeftest(mod15_8FE, vcov = vcovSCC(mod15_8FE, type = "HC1"))


# =============================================================================
#' # Part 7. The random effects (RE) estimator
# =============================================================================
#' The random effects estimator treats the individual differences as being randomly assigned to
#' the individuals. Rather than estimate them as parameters as we did in the fixed effects model,
#' here they are incorporated into the model's error, which in a panel will have a specific structure.
#' One of the key advantages of the random effects model is that parameters
#' on time invariant regressors can be estimated.
#' Key assumption: the unobserved individual effects u_i are uncorrelated with the regressors.
#' The parameter estimates are actually obtained through feasible generalized least squares.
#' The transformation that is used on the variables of the model is sometimes referred to as quasi-demeaning.
#' The transformation parameter is called theta in the plm output.
#' As theta -> 1, RE approaches FE; as theta -> 0, RE approaches pooled OLS.

# -----------------------------------------------------------------------------
#' ## Example 15.9 (POE5)
#' For the 1000 Chinese chemical firms the random effects are estimated as:
mod15_9RE <- plm(lsales ~ lcapital + llabor, data = chem3PDF, model = "random")
summary(mod15_9RE)
#' The estimated value of alpha (the quasi-demeaning weight) in the random effects
#' estimator is 0.7353. plm calls this theta.


# =============================================================================
#' # Part 8. The Hausman test (FE vs. RE)
# =============================================================================
#' If the random individual effects are correlated with regressors,
#' then the random effects estimator will not be consistent.
#' A statistical test of this proposition should be done whenever this estimator
#' is used in order to reduce the chance of model misspecification.
#'
#' H0: The RE model is consistent (u_i uncorrelated with the regressors).
#' H1: The RE model is inconsistent.

# -----------------------------------------------------------------------------
#' ## Hausman test for endogeneity, Example 15.12 (POE5)
phtest(mod15_8FE, mod15_9RE)
#' Decision rule:
#' - Reject H0: RE is inconsistent, use the Fixed Effects (FE) model.
#' - Do not reject H0: RE is consistent and more efficient, so RE can be preferred
#'   (gives precision and allows time-invariant regressors).


# =============================================================================
#' # Part 9. A complete workflow: the wage equation (nls_panel)
# =============================================================================
#' We now go through all the estimators on the wage data from Example 15.1,
#' and add some diagnostic tests.

# -----------------------------------------------------------------------------
#' ## A pooled wage model
#' Note: This model is not considered in the POE5 text (see POE4 Chapter 15).
#' Coefficients are partial associations under the strong assumption of no omitted
#' individual heterogeneity. Always use at least cluster-robust SEs at the id level.
wage.pooled <- plm(lwage ~ educ + exper + I(exper^2) + tenure + I(tenure^2) + black + south + union,
                   model = "pooling", data = nlsPDF)

pooled.mod <- wage.pooled %>% tidy(.) # output
pooled.mod
pooled.mod.robust <- tidy(coeftest(wage.pooled,
                                   vcov = vcovHC(wage.pooled, type = "HC0", cluster = "group")))

#' Conventional and cluster-robust results side by side
pooled.mod %>% left_join(select(pooled.mod.robust, -estimate), by = "term") %>%
  kable(., align = "c", digits = 3,
        col.names = c("Variable", "Coefficient",
                      "Std. Error OLS", "t-value OLS", "p-value OLS",
                      "Std. Error robust", "t-value robust", "p-value robust"))

#' Note that a plain OLS model would give you the same estimates as plm(model = "pooling").
wage.ols <- lm(lwage ~ educ + exper + I(exper^2) + tenure + I(tenure^2) + black + south + union,
               data = nls_panel)
pooled.mod.ols <- wage.ols %>% tidy(.)
cbind(pooled.mod$estimate, pooled.mod.ols$estimate)
all.equal(pooled.mod$estimate, pooled.mod.ols$estimate)

# -----------------------------------------------------------------------------
#' ## The fixed effects model: LSDV with a small N (POE4)
#' With only N=10 women we can look at all the dummy variable coefficients.
nls10 <- pdata.frame(nls_panel[nls_panel$id %in% 1:10, ], index = c("id", "year"))
pdim(nls10)

#' Dummy variable estimator for small N, using OLS
wage.fixed.10.ols <- lm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + union + factor(id) - 1,
                        data = nls10)
kable(tidy(wage.fixed.10.ols), digits = 4)

#' Alternative restricted model, with only one intercept, POE4 Table 15.4
wage.fixed.pooled.10 <- plm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + union,
                            data = nls10, model = "pooling")
kable(tidy(wage.fixed.pooled.10), digits = 4)

#' It is not necessary to use the OLS dummy variable approach,
#' use the option model = "within" in plm(). POE4 Table 15.6
wage.fixed.10 <- plm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + union,
                     data = nls10, model = "within")
kable(tidy(wage.fixed.10), digits = 4)

#' Recover the 10 intercepts; compare with the factor(id) coefficients above
fixef(wage.fixed.10, type = "level")

#' Test of poolability (POE4 p. 546)
pFtest(wage.fixed.10, wage.fixed.pooled.10)

# -----------------------------------------------------------------------------
#' ## The fixed effects model on the complete panel, Table 15.7 (POE5)
#' Time-invariant variables (educ, black) cannot be included in the FE model.
wage.within <- plm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + south + union,
                   data = nlsPDF, model = "within")
wage.within %>% tidy() %>% kable(., digits = 5)

#' The LSDV estimator gives the same slopes (716 dummies, so we print the slopes only)
wage.fixed.ols <- lm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + south + union + factor(id) - 1,
                     data = nls_panel)
tidy(wage.fixed.ols) %>%
  filter(term %in% c("exper", "I(exper^2)", "tenure", "I(tenure^2)", "south", "union")) %>%
  kable(., digits = 5)

#' Histogram of all 716 intercepts
wage.within %>% fixef() %>% histogram(., width = 0.1)

#' Test of poolability: FE vs. pooled OLS with the same regressors
wage.fixed.pooled <- plm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + south + union,
                         data = nlsPDF, model = "pooling")
pFtest(wage.within, wage.fixed.pooled)
#' Rejecting H0 implies worker-specific intercepts matter, so FE is preferred to pooled OLS.

# -----------------------------------------------------------------------------
#' ## The random effects model
#' Breusch-Pagan LM test for random effects (POE4 p. 554).
#' H0: no individual effects (the variance of u_i is zero), i.e., pooled OLS is fine.
plmtest(wage.pooled, effect = "individual", type = "bp")
#' Rejecting H0 means pooled OLS is inappropriate; use FE or RE.

#' Now we can also include the time invariant variables
wage.random <- plm(lwage ~ educ + exper + I(exper^2) + tenure + I(tenure^2) + black + south + union,
                   data = nlsPDF, random.method = "swar", model = "random")
kable(tidy(wage.random), digits = 4)
summary(wage.random)  # shows theta

#' Hausman test: FE vs. RE (compares the common coefficients)
phtest(wage.within, wage.random)
#' The null hypothesis, that the individual random effects are exogenous, is rejected.
#' This makes the random effects estimator inconsistent.
#' The fixed effects model is the preferred model.

# -----------------------------------------------------------------------------
#' ## Further panel diagnostics
#' These help students see why we choose certain standard errors or include time effects.

#' 1) Serial correlation in FE residuals (Wooldridge test)
#' H0: no serial correlation. Reject: cluster by group and/or use Driscoll-Kraay SEs.
pwartest(wage.within)

#' 2) Cross-sectional dependence (Pesaran CD test)
#' H0: cross-sectional independence. Reject: consider Driscoll-Kraay (vcovSCC).
pcdtest(wage.within, test = "cd")

#' 3) Time fixed effects: do we need to control for common shocks in each year?
#' H0: time effects are jointly zero. If rejected, keep factor(year).
wage.within.tw <- plm(lwage ~ exper + I(exper^2) + tenure + I(tenure^2) + south + union + factor(year),
                      data = nlsPDF, model = "within")
pFtest(wage.within.tw, wage.within)
#' The same can be done with model = "within", effect = "twoways" (individual and time effects).


# =============================================================================
#' # Part 10. Within and between effects
# =============================================================================

# -----------------------------------------------------------------------------
#' ## Example 15.14 but as a random effects within-between model
#' The within-between model decomposes each time-varying regressor into
#' a within part (x_it - xbar_i) and a between part (xbar_i).
#' - Within coefficient: short-run (intra-firm) elasticity. Equal to the FE estimate.
#' - Between coefficient: cross-sectional (inter-firm) elasticity.
#' If they differ materially, short-run and long-run responses differ, and the
#' RE assumption (u_i uncorrelated with the regressors) is doubtful.
# browseURL("https://panelr.jacob-long.com/articles/wbm")

chem3PanelR <- panel_data(chemical3, id = firm, wave = year)
names(chem3PanelR)

#' Syntax: dependent ~ time_varying | time_invariant | cross_level_interactions
#' The time-varying variables go in the first part; wbm() creates the within and
#' between parts automatically. We have no time-invariant variables here.
model.wb <- wbm(lsales ~ lcapital + llabor, data = chem3PanelR)
summary(model.wb)
#' Compare the within effects with the FE estimates in mod15_8FE:
coef(mod15_8FE)

# -----------------------------------------------------------------------------
#' ## Plots of within vs. between variation (wage data)
#' These plots build intuition for what the FE and RE models are actually doing.

#' Prepare data for plotting by demeaning and averaging (drop the pseries class first)
plot_df <- as.data.frame(nlsPDF) %>%
  mutate(id_num = as.integer(as.character(id)),
         exper_within = exper - ave(exper, id_num))   # deviation from person-mean

nls_means <- plot_df %>%
  group_by(id_num) %>%
  summarise(lwage_mean = mean(lwage), exper_mean = mean(exper), .groups = "drop")

#' 1) Within plot: how wages change as an individual's experience deviates from her mean.
#' The slope corresponds to the fixed effects idea.
ggplot(plot_df, aes(x = exper_within, y = lwage, group = id_num)) +
  geom_line(alpha = 0.10, linewidth = 0.2) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE) +
  labs(title = "Within effect: deviations from each worker's mean experience",
       x = "Experience (demeaned by worker)", y = "Log wage")

#' 2) Between plot: cross-sectional relationship between workers' averages.
ggplot(nls_means, aes(x = exper_mean, y = lwage_mean)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Between effect: averages across workers",
       x = "Average experience", y = "Average log wage")

#' Interpretation: compare the slopes. If the within slope differs from the between slope,
#' the return to one more year of experience for one person (within) differs from the
#' wage gap between two people who differ by one year in average experience (between).


# =============================================================================
#' # Part 11. The Hausman-Taylor estimator
# =============================================================================
#' FE is consistent but cannot estimate the effect of time-invariant variables (like educ).
#' RE can, but is inconsistent if some regressors are correlated with u_i.
#' The Hausman-Taylor estimator is an instrumental variables (IV) version of RE:
#' it uses the exogenous variables in the model as instruments for the endogenous ones.
#'
#' Formula syntax in plm (three parts):
#'   y ~ all regressors | all exogenous regressors | time-varying endogenous regressors
#' Here we treat educ (time-invariant) and south (time-varying) as endogenous.
#' Note: older code used model = "ht" or pht(). These are deprecated in recent plm versions.
wage.HT <- plm(lwage ~ educ + exper + I(exper^2) + tenure + I(tenure^2) + black + south + union |
                 exper + I(exper^2) + tenure + I(tenure^2) + union + black |
                 south,
               data = nlsPDF, model = "random", random.method = "ht", inst.method = "baltagi")
summary(wage.HT)
kable(tidy(wage.HT), digits = 4)

#' Compare the return to education in RE and Hausman-Taylor
coef(wage.random)["educ"]
coef(wage.HT)["educ"]


# =============================================================================
#' # Part 12. Sets of regression equations (SUR), POE4 Chapter 15
# =============================================================================
#' When T is large relative to N, we can allow every individual to have its own
#' intercept AND its own slopes, i.e., estimate one equation per individual.
#' Investment data for General Electric and Westinghouse (grunfeld2 in POE4).
#' We use the Grunfeld data from the AER package, keeping the same two firms.
#' The numbers may differ slightly from the POE4 textbook data.
data("Grunfeld", package = "AER")
grunfeld2 <- subset(Grunfeld, firm %in% c("General Electric", "Westinghouse"))
grunfeld2$firm <- factor(grunfeld2$firm, levels = c("General Electric", "Westinghouse"),
                         labels = c("GE", "WE"))
PDgrun <- pdata.frame(grunfeld2, index = c("firm", "year"))
pdim(PDgrun)
kable(head(PDgrun), align = "c")

#' Pooled model: same intercept and slopes for both firms
grun.pool <- plm(invest ~ value + capital, model = "pooling", data = PDgrun)
kable(tidy(grun.pool), digits = 5)

#' Firm-specific intercepts and slopes (interactions with firm)
grun.fe <- plm(invest ~ value*firm + capital*firm, model = "pooling", data = PDgrun)
kable(tidy(grun.fe), digits = 5)

#' F-test: are intercepts and slopes equal across firms?
pFtest(grun.fe, grun.pool)

#' Two separate regressions
grun1.pool <- lm(invest ~ value + capital, data = grunfeld2, subset = firm == "GE")
grun2.pool <- lm(invest ~ value + capital, data = grunfeld2, subset = firm == "WE")
glance(grun1.pool)
glance(grun2.pool)

#' Goldfeld-Quandt test. The data are sorted by firm, so the sample is split
#' between GE (first half) and WE (second half).
#' Null hypothesis: the error variances are equal.
gqtest(invest ~ value + capital, point = 0.5, alternative = "two.sided", data = grunfeld2)
#' Rejected: the separate regressions have different error variances.

#' Seemingly Unrelated Regressions (SUR) with systemfit.
#' Note that the data have to be a pdata.frame. systemfit then estimates one
#' equation per firm and allows the errors of the two equations to be correlated
#' (contemporaneous correlation).
grunf.SUR <- systemfit(invest ~ value + capital, method = "SUR", data = PDgrun)
summary(grunf.SUR)

#' The coefficient names are needed for the hypothesis test below
coef(grunf.SUR)

#' Testing cross-equation hypothesis: same intercept and slopes in both equations
RMatrix <- c("GE_(Intercept) = WE_(Intercept)", "GE_value = WE_value", "GE_capital = WE_capital")
linearHypothesis(grunf.SUR, RMatrix)

#' Same coefficients, pooled model
grunf.SUR.pooled <- systemfit(invest ~ value + capital, method = "SUR",
                              data = PDgrun, pooled = TRUE)
summary(grunf.SUR.pooled)


# =============================================================================
#' # Part 13. Applied exercises
# =============================================================================

# -----------------------------------------------------------------------------
#' ## POE5 Exercise 15.18 (POE4 Exercise 15.6): Mexican sex worker data
#' Data definition file
# browseURL("http://www.principlesofeconometrics.com/poe5/data/def/mexican.def")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mexican.rdata"))

str(mexican)
mexicanPDF <- pdata.frame(mexican, index = c("id", "trans"))
pdim(mexicanPDF)
names(mexicanPDF)

#' 1) sex worker characteristics: i) "age", ii) "attractive" iii) "school"
#' 2) client characteristics: i) "regular", ii) "rich", iii) "alcohol"
#' 3) transaction characteristics: i) "lnprice", ii) "nocondom", iii) "bar" & "street" (basis is "othersite")
mexican %>% summarise(mean(bar))
mexican %>% summarise(mean(street))
mexican %>% summarise(mean(othersite))

#' lnprice    = logarithm of price of transaction
#' bar        = 1 if transaction originated in a bar; 0 otherwise
#' street     = 1 if transaction originated in a street; 0 otherwise
#' othersite  = 1 if transaction originated in another site; 0 otherwise
#' nocondom   = 1 if a condom was not used; 0 otherwise
#' attractive = 1 if the sex worker is attractive; 0 otherwise
#' school     = 1 if sex worker has completed secondary school or higher; 0 otherwise
#' age        = age of sex worker in years
#' rich       = 1 if client is rich; 0 otherwise
#' regular    = 1 if client is regular; 0 otherwise
#' alcohol    = 1 if client consumed alcohol prior to the transaction

#' Abstract:
#' While condoms are an effective defense against the transmission of
#' HIV, large numbers of sex workers are not using them. We argue that
#' some sex workers are willing to take the risk because clients are willing
#' to pay more to avoid using condoms. Using data from Mexico, we
#' estimate that sex workers received a 23 percent premium for unprotected sex.

#' a. OLS model
lnprice.ols <- lm(lnprice ~ age + attractive + school + regular + rich + alcohol + nocondom + bar + street,
                  data = mexican)
summary(lnprice.ols)

#' Age: log-linear model, 1 unit increase in age (1 year)
coef(lnprice.ols)["age"]*100 # percent lower price

#' Dummy variables in a log-linear model: the percentage effect is 100*(exp(b)-1).
#' Attractive
100*(exp(coef(lnprice.ols)["attractive"]) - 1) # percent higher price
#' School
100*(exp(coef(lnprice.ols)["school"]) - 1)
#' Regular client
100*(exp(coef(lnprice.ols)["regular"]) - 1)
#' Rich client
100*(exp(coef(lnprice.ols)["rich"]) - 1)
#' Alcohol
100*(exp(coef(lnprice.ols)["alcohol"]) - 1)
#' Nocondom
100*(exp(coef(lnprice.ols)["nocondom"]) - 1)
#' Note: the "premium" for not using a condom is small and statistically insignificant.
#' Bar
100*(exp(coef(lnprice.ols)["bar"]) - 1)
#' Street
100*(exp(coef(lnprice.ols)["street"]) - 1)

#' All the coefficients are significantly different from zero at the
#' 5% level except that of NOCONDOM. The signs are plausible.

#' b. The 95% interval estimate of the coefficient of nocondom is
ci <- confint(lnprice.ols)
ci["nocondom", ]
#' That is, we estimate that the risk premium is between
100*(exp(ci["nocondom", ]) - 1)
#' The interval covers zero, which means that the coefficient is not significantly
#' different from zero at the 5% level. We conclude that we have not estimated
#' the risk premium precisely, if there is one.

#' c. The unobserved characteristics of the sex worker that may be correlated with the explanatory
#' variables might include how much income the sex worker has other than from sex work, the
#' amount of experience the worker has, the type of services offered, whether or not she has a
#' sexually transmitted disease, whether she has a partner and/or has children, whether she smokes
#' or uses illegal drugs or uses alcohol.

#' d. Fixed Effects Model
#' The sex worker characteristics are time invariant, meaning that the fixed effects,
#' the individual specific indicator variables, are perfectly correlated with them.
#' Or, equivalently, the within-transformed values of these variables are zero.
#' Their effect cannot be separated from the individual effects.
lnprice.fixed <- plm(lnprice ~ regular + rich + alcohol + nocondom + bar + street,
                     data = mexicanPDF, model = "within")
summary(lnprice.fixed)
#' All of the coefficients in the fixed effects model are statistically significant at the 5% level.

stargazer(lnprice.ols, lnprice.fixed, type = "text")

#' Intercepts
fixef(lnprice.fixed, type = "level") %>% histogram(., width = 0.1)
summary(fixef(lnprice.fixed, type = "level"))
exp(3)
exp(8)

#' Mean of all individual intercepts, and in levels
mean(fixef(lnprice.fixed, type = "level"))
exp(mean(fixef(lnprice.fixed, type = "level")))

#' e. Do the F-test (equation 15.20)
glance(lnprice.ols)   # SSER
glance(lnprice.fixed) # SSEU

SSER <- sum(resid(lnprice.ols)^2)
SSEU <- sum(resid(lnprice.fixed)^2)
N  <- length(fixef(lnprice.fixed))
NT <- nrow(mexicanPDF)
K  <- length(coef(lnprice.fixed))
#' F-value
((SSER - SSEU)/(N - 1)) / (SSEU/(NT - N - K))
#' Critical F
qf(0.99, N - 1, NT - N - K)
#' We reject the null hypothesis that there are no individual differences among sex workers.

#' f. Percentage effects in the FE model
100*(exp(coef(lnprice.fixed)) - 1)
#' Note: the "premium" for not using a condom is now statistically significant.

#' g. The 95% interval estimate of the coefficient of nocondom is
ci <- confint(lnprice.fixed)
ci["nocondom", ]
#' That is, we estimate that the risk premium is between
100*(exp(ci["nocondom", ]) - 1)
#' The interval does not cover zero, which means that the coefficient is significantly
#' different from zero at the 5% level. This interval is slightly narrower than the interval based
#' on the OLS estimates, and of course only covers positive values, suggesting that there is a risk
#' premium paid for unprotected sex.

#' --------- End of the exercise, extensions below ------------------------

#' Random Effects Model
lnprice.random <- plm(lnprice ~ age + attractive + school + regular + rich + alcohol + nocondom + bar + street,
                      data = mexicanPDF, random.method = "swar", model = "random")
summary(lnprice.random)

stargazer(lnprice.ols, lnprice.fixed, lnprice.random, type = "text")

#' Treating the effects as random instead of fixed and adding the sex worker characteristics
#' has had a dramatic effect on some of the common coefficients.
#' We pick the common coefficients by name:
common <- names(coef(lnprice.fixed))
round(100*(exp(coef(lnprice.fixed)) - 1), 1)
round(100*(exp(coef(lnprice.random)[common]) - 1), 1)

#' Difference in effect between the two models
round(100*(exp(coef(lnprice.fixed)) - 1), 1) - round(100*(exp(coef(lnprice.random)[common]) - 1), 1)

#' Price premium on unprotected sex with an educated attractive woman (RE model)
round(100*(exp(coef(lnprice.random)[c("nocondom", "attractive", "school")]) - 1), 1)
round(100*(exp(sum(coef(lnprice.random)[c("nocondom", "attractive", "school")])) - 1), 1)

#' Hausman test (POE4 section 15.5.3), all common coefficients jointly
phtest(lnprice.fixed, lnprice.random)

#' Hausman test on each coefficient separately
bFE <- tidy(lnprice.fixed) %>% select(term, estimate, std.error)
bRE <- tidy(lnprice.random) %>% select(term, estimate, std.error) %>% filter(term %in% common)
df <- merge(bFE, bRE, by = "term")

#' Results for the Hausman test on each difference between the fixed effects and random
#' effects estimates are given in the following table. The test uses the large sample
#' normal distribution.
df %>% mutate(b  = estimate.x - estimate.y,
              se = sqrt(std.error.x^2 - std.error.y^2),
              t  = b/se,
              p  = round(2*pnorm(-abs(t)), 4)) %>%
  select(term, b, se, t, p)
#' At a 5% level of significance, there is a significant difference between all coefficients
#' except those for BAR. Thus, we reject a null hypothesis that the individual random effects
#' are uncorrelated with the variables in the model. The fixed effects estimates are more
#' reliable in this instance because they are consistent.

#' The Hausman-Taylor estimator
#' If a sex worker has individual characteristics that make her a risk taker, or, conversely,
#' risk averse, then NOCONDOM is likely to be correlated with the individual effect.
#' We treat NOCONDOM as endogenous; all other variables are exogenous.
lnprice.HT <- plm(lnprice ~ bar + street + nocondom + rich + regular + alcohol + attractive + school + age |
                    bar + street + rich + regular + alcohol + attractive + school + age |
                    nocondom,
                  data = mexicanPDF, model = "random", random.method = "ht", inst.method = "baltagi")
summary(lnprice.HT)

#' The results are very similar to those obtained with FE. There have been no dramatic
#' changes in the coefficient estimates and REGULAR, ALCOHOL and STREET continue to be
#' insignificant at a 5% level of significance.
#' In this case, the extra percentage premium for having unprotected sex with an attractive
#' secondary-educated sex worker, compared with protected sex with an unattractive
#' uneducated sex worker is:
bHT <- coef(lnprice.HT)[c("nocondom", "attractive", "school")]
bHT
100*(exp(sum(bHT)) - 1)

# -----------------------------------------------------------------------------
#' ## POE5 Exercise 15.22: Crime in North Carolina counties
#' Data definition file
# browseURL("http://www.principlesofeconometrics.com/poe5/data/def/crime.def")
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/crime.rdata"))

str(crime)
crimePDF <- pdata.frame(crime, index = c("county", "year"))
pdim(crimePDF)

#' a. **A priori assumptions:**
#' (i) If deterrence increases crime rates should drop.
#' (ii) If wages in the private sector increase the return to legal activities increases relative to the
#' return to illegal activities. Therefore crime rates should drop.
#' (iii) Higher population density should be linked with a higher residential crime rate.
#' (iv) Young males are the most likely demographic group to be involved in illegal activities.
#' Thus, an increase in the percentage of young males should increase the crime rate.

#' b. **A pooled OLS model**
#' lcrmrte   log(crimes committed per person)
#' lprbarr   log('probability' of arrest)
#' lprbconv  log('probability' of conviction)
#' lprbpris  log('probability' of prison sentence)
#' lavgsen   log(avg. sentence, days)
#' lwmfg     log(weekly wage, manufacturing)
noquote(names(crime))
lcrmrte.ols <- lm(lcrmrte ~ lprbarr + lprbconv + lprbpris + lavgsen + lwmfg, data = crime)
summary(lcrmrte.ols)

#' This is a log-log model, hence all parameters are interpreted as elasticities.
#' Variables measuring deterrence of the legal system:
#' lprbarr : when the 'probability' of arrest increases by 1% the crimes committed per person (crime rate)
#' changes by `r round(coef(lcrmrte.ols)[2],2)` %.
#' lprbconv : when the 'probability' of conviction increases by 1% the crime rate
#' changes by `r round(coef(lcrmrte.ols)[3],2)` %.
#' lprbpris : when the 'probability' of prison sentence increases by 1% the crime rate
#' changes by `r round(coef(lcrmrte.ols)[4],2)` %. This is somewhat surprising.
#' lavgsen : when the avg. sentence in days increases by 1% the crime rate
#' changes by `r round(coef(lcrmrte.ols)[5],2)` %.
#' If wages, measured by lwmfg (weekly wage, manufacturing), increase by 1% the crime rate
#' changes by `r round(coef(lcrmrte.ols)[6],2)` %. This is somewhat surprising.
#'
#' For the variables that describe the deterrence effect of the legal system we would expect that the
#' coefficients would be negative. We find that all of these coefficients are negative except for the
#' coefficient of LPRBPRIS. The variable LWMFG, which represents wages in the private sector, has a
#' positive coefficient that is not consistent with our expectations.
#' All coefficients are significantly different from zero at a 5% level of significance except for the
#' coefficient of LAVGSEN.

#' c. **The Fixed Effects model**
lcrmrte.fixed <- plm(lcrmrte ~ lprbarr + lprbconv + lprbpris + lavgsen + lwmfg, data = crimePDF, model = "within")
summary(lcrmrte.fixed)

#' Both models in the same output:
stargazer(lcrmrte.ols, lcrmrte.fixed, type = "text")

#' All estimated coefficients have the expected sign except for LAVGSEN. Moreover, all
#' estimated coefficients are significantly different from zero at a 5% level of significance except
#' for the coefficient for LAVGSEN (the avg. sentence, days).
#'
#' (ii) The coefficient on LPRBARR suggests that a 1% increase in the probability of being
#' arrested results in a `r round(coef(lcrmrte.fixed)[1],2)` % change in the crime rate.
#' This estimated elasticity is about one third (`r round(coef(lcrmrte.fixed)[1]/coef(lcrmrte.ols)[2],2)`) of the
#' estimated elasticity from the pooled OLS model. Thus, once we allow for county heterogeneity, the
#' deterrent effect of being arrested is much less.
#'
#' (iii) The coefficient on LAVGSEN suggests that a 1% increase in the average prison sentence
#' results in a `r round(coef(lcrmrte.fixed)[4],2)` % change in the crime rate.
#' However, a two tail t-test on the significance of this estimate yields a t-statistic of
#' `r round(broom::tidy(lcrmrte.fixed)[4,4],3)` and a p-value of `r round(broom::tidy(lcrmrte.fixed)[4,5],3)`.
#' Thus, a null hypothesis that the coefficient of LAVGSEN is zero is not rejected. There is no
#' support for the idea that longer prison sentences are a deterrent to the crime rate.
#' (The `r ...` expressions are filled in when the script is rendered with knitr::spin.)

#' d. **F-test on equal intercepts in the Fixed Effects model**
#' This is a histogram of all the 90 intercepts:
fixef(lcrmrte.fixed, type = "level") %>% histogram(., width = 0.1)

#' Do the F-test (equation 15.20)
SSER <- sum(resid(lcrmrte.ols)^2)
SSEU <- sum(resid(lcrmrte.fixed)^2)
N  <- length(fixef(lcrmrte.fixed))
NT <- nrow(crimePDF)
K  <- length(coef(lcrmrte.fixed))
#' F-value
((SSER - SSEU)/(N - 1)) / (SSEU/(NT - N - K))
#' Critical F
qf(0.95, N - 1, NT - N - K)
#' p-value
pf(((SSER - SSEU)/(N - 1)) / (SSEU/(NT - N - K)), df1 = N - 1, df2 = NT - N - K, lower.tail = FALSE)

#' To test H0: beta_{1,1} = beta_{1,2} = ... = beta_{1,90} against the alternative that not all of the
#' intercepts are equal, we use the usual F-test for testing a set of linear restrictions.
#' We reject H0 and conclude that the county level effects are not all equal.

#' e. **Policy implications**
#' According to the fixed effects estimates, the explanatory variables which have the expected
#' signs and a significant effect on the crime rate are LPRBARR (the probability of being
#' arrested), LPRBCONV (the probability of being convicted), LPRBPRIS (the probability of prison sentence)
#' and LWMFG (weekly wage, manufacturing). Out of these variables, those that have the largest effect
#' on the crime rate, and are reasonable to implement as public policy, will be the most effective in
#' dealing with crime. Improving policing and court policies that increase the probability of arrest,
#' conviction and imprisonment are likely to be effective, but lengthening the term of
#' imprisonment is not (LAVGSEN). Opportunities for higher wages are also likely to be a productive
#' direction for public policy.

#' **Why are** ldensity (log people per sq. mile) **and** lpctymle (log percent young male)
#' **not included in the model, when they are mentioned in a.?**
ols.mod2   <- update(lcrmrte.ols,   . ~ . + ldensity + lpctymle)
fixed.mod2 <- update(lcrmrte.fixed, . ~ . + ldensity + lpctymle)

#' Both updated models in the same output:
stargazer(ols.mod2, fixed.mod2, type = "text")

# -----------------------------------------------------------------------------
#' ## POE5 Exercise 15.29: Panel IV with one endogenous variable (police per capita)
#' Police per capita (lpolpc) is likely endogenous: more crime leads to more police.
#' Instruments: ltaxpc (tax revenue per capita) and lmix (offense mix).

#' a. Estimate the reduced form (1st stage) of the endogenous variable, with the two instruments
lpolpc.red_form <- lm(lpolpc ~ lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                      data = crime)
summary(lpolpc.red_form)

#' Test of joint significance of the two IV
linearHypothesis(lpolpc.red_form, c("ltaxpc=0", "lmix=0"))
tidy(lpolpc.red_form)

#' The two instruments have positive coefficients and are significant individually.
#' The F-test of their joint significance is
#' `r linearHypothesis(lpolpc.red_form, c("ltaxpc=0", "lmix=0"))$F[2]` > 10, the rule of thumb value.
#' Thus we can reject the null hypothesis that the IV are weak using this criterion.

#' b. 2SLS Estimation
lcrmrte.2sls <- tsls(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban,
                     ~ lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                     data = crime)
summary(lcrmrte.2sls)

#' The deterrence variables, the log of the probability of arrest (LPRBARR),
#' the log of the probability of conviction (LPRBCONV),
#' the log of average prison sentence (LAVGSEN), all have negative
#' and significant coefficients, indicating that they are having the desired effect.
#' The log of police per capita has a positive coefficient
#' and is significant at the 10% level.

#' c. Testing for exogeneity
#' http://eclr.humanities.manchester.ac.uk/index.php/IV_in_R
#' Regression based Hausman test: add the first stage residuals to the OLS model
lcrmrte.ols_w_1stage_res <- lm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban +
                                 resid(lpolpc.red_form), data = crime)
summary(lcrmrte.ols_w_1stage_res)

#' Hausman-Wu test
HausWutest <- waldtest(lcrmrte.ols_w_1stage_res, . ~ . - resid(lpolpc.red_form))
print(HausWutest)

#' Including the first stage residuals into the equation and estimating
#' it by OLS we find that the first stage residuals have a coefficient of
#' 0.1048 and a t = 0.78 (F = 0.6083 with p = 0.4358).
#' Thus based on this test we do not find evidence that LPOLPC is endogenous.

#' The Sargan test of the validity of the surplus instrument yields:
Sargan_reg <- lm(resid(lcrmrte.2sls) ~ lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                 data = crime)
Sargan_reg_sm <- summary(Sargan_reg)

Sargan_test <- Sargan_reg_sm$r.squared*nrow(crime)
print(Sargan_test)
print(1 - pchisq(Sargan_test, 1))  # prints p-value
#' Thus we reject the validity of the surplus IV, making the 2SLS results suspect at best.

#' d. The Reduced Form Fixed Effects model
lpolpc.red_form.fix <- plm(lpolpc ~ lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                           data = crimePDF, model = "within")
summary(lpolpc.red_form.fix)

#' Both of the instruments are significant at the 5% level.
#' Their joint test of significance yields an F:
linearHypothesis(lpolpc.red_form.fix, c("ltaxpc=0", "lmix=0"))
#' which indicates they are statistically significant.

#' e. Fixed Effects 2SLS. In plm, the instruments come after the "|".
lcrmrte.2sls.fix <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban
                        | lprbarr + lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                        data = crimePDF, model = "within")
summary(lcrmrte.2sls.fix)
#' The deterrent variables behave as before, with the log of per capita police still
#' positive and significant.

#' f. Hausman test in the FE model
lcrmrte.fixed_w_1stage_res <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban +
                                    resid(lpolpc.red_form.fix),
                                  data = crimePDF, model = "within")
summary(lcrmrte.fixed_w_1stage_res)

HausWutest.fixed <- waldtest(lcrmrte.fixed_w_1stage_res, . ~ . - resid(lpolpc.red_form.fix))
print(HausWutest.fixed)

#' Thus using the within data we find (weak) endogeneity of the log of police
#' per capita and we do not reject the validity of the surplus IV.
#' It should be noted that these tests were not designed for use with panel data.

# -----------------------------------------------------------------------------
#' ## POE5 Exercise 15.30: Panel IV with two endogenous variables
#' Now both police per capita (lpolpc) and the probability of arrest (lprbarr)
#' are treated as endogenous.
lpolpc.red_form <- lm(lpolpc ~ lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                      data = crime)
summary(lpolpc.red_form)

lprbarr.red_form <- lm(lprbarr ~ lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                       data = crime)
summary(lprbarr.red_form)

#' Test of joint significance of the two IV
linearHypothesis(lpolpc.red_form,  c("ltaxpc=0", "lmix=0"))
linearHypothesis(lprbarr.red_form, c("ltaxpc=0", "lmix=0"))

#' Note that for LPOLPC both instruments are significant and have a joint
#' F-statistic value of 20.835. For LPRBARR only LMIX is significant but the
#' joint F-test of significance is 113.52.
#' Recall from page 504 of POE5 that when there are 2
#' endogenous variables the two F-tests are not adequate to measure IV strength.
#' An alternative is discussed in Appendix 10A, the minimum eigenvalue statistic
#' or the Cragg-Donald F-statistic.
#' We can reject the null hypothesis that the IV are weak.

#' Hausman test
lcrmrte.ols_w_1stage_res <- lm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban +
                                 resid(lpolpc.red_form) + resid(lprbarr.red_form), data = crime)
summary(lcrmrte.ols_w_1stage_res)

#' The regression based on the Hausman test shows that the first stage
#' residuals are significant. The F-test of their joint significance is:
linearHypothesis(lcrmrte.ols_w_1stage_res, c("resid(lpolpc.red_form)=0", "resid(lprbarr.red_form)=0"))
#' We conclude that one or both of the variables LPOLPC and LPRBARR is endogenous.

#' 2SLS Estimation
lcrmrte.2sls <- tsls(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg + west + urban,
                     ~ lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                     data = crime)
summary(lcrmrte.2sls)
#' We see that LPOLPC and LPRBARR are insignificant.
#' LPRBCONV and LAVGSEN have negative coefficients and the
#' coefficient of the probability of conviction variable is significant.

#' Fixed Effects reduced forms
lpolpc.fixed.red_form <- plm(lpolpc ~ lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                             data = crimePDF, model = "within")
summary(lpolpc.fixed.red_form)

lprbarr.fixed.red_form <- plm(lprbarr ~ lprbconv + lavgsen + lwmfg + ltaxpc + lmix,
                              data = crimePDF, model = "within")
summary(lprbarr.fixed.red_form)

#' We see that LMIX is a strong IV in both cases. The tests of joint significance are:
linearHypothesis(lpolpc.fixed.red_form,  c("ltaxpc=0", "lmix=0"), test = "F")
linearHypothesis(lprbarr.fixed.red_form, c("ltaxpc=0", "lmix=0"), test = "F")
#' We reject the null hypothesis that the IV are weak.

#' Fixed Effects 2SLS
lcrmrte.fixed.2sls <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg
                          | lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix,
                          data = crimePDF, model = "within")
summary(lcrmrte.fixed.2sls)

#' e. The results are once again not what we would expect, as LPOLPC has a positive
#' and significant coefficient. LAVGSEN is now insignificant and the other
#' two deterrence variables have negative coefficients.

#' f. Hausman test in the FE-2SLS model
lcrmrte.fixed.2sls_w_1st <- plm(lcrmrte ~ lpolpc + lprbarr + lprbconv + lavgsen + lwmfg +
                                  resid(lpolpc.fixed.red_form) + resid(lprbarr.fixed.red_form)
                                | lprbconv + lavgsen + lwmfg + west + urban + ltaxpc + lmix +
                                  resid(lpolpc.fixed.red_form) + resid(lprbarr.fixed.red_form),
                                data = crimePDF, model = "within")
summary(lcrmrte.fixed.2sls_w_1st)
linearHypothesis(lcrmrte.fixed.2sls_w_1st, c("resid(lpolpc.fixed.red_form)=0",
                                             "resid(lprbarr.fixed.red_form)=0"))
#' We find slight evidence of endogeneity.

#' ---------------------------- End of script ---------------------------------

