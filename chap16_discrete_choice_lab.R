#' ===========================================================================
#' # Chapter 16: Qualitative and Limited Dependent Variable Models
#' ## R coding session -- a standalone resource
#' ===========================================================================
#'
#' Main reference: Hill, Griffiths & Lim, Principles of Econometrics, 5th ed.
#' (POE5), Chapter 16. Page numbers refer to POE5 unless marked POE4.
#' Some marketing examples follow Franses & Paap, Quantitative Models in
#' Marketing Research.
#'
#' The script moves from simple to more complicated models:
#'
#'   Part 1  Binary choice: LPM, probit and logit (transport data)
#'           - The linear probability model
#'           - Example 16.3: probit by maximum likelihood "by hand"
#'           - Probit: predicted probabilities, marginal effects, elasticities
#'           - The logistic function and logit elasticities
#'           - Logit and comparison with probit
#'   Part 2  Binary choice with several regressors (coke data, Example 16.6)
#'           - Hypothesis tests (Wald, Example 16.7)
#'           - Predicted probabilities, marginal effects, elasticities
#'   Part 3  Binary choice in marketing (Heinz vs Hunts, Franses & Paap)
#'           - Model evaluation: pseudo R2, confusion matrix, ROC, hold-out
#'   Part 4  Interaction terms in a logit model (Cracker data)
#'   Part 5  Multinomial logit (MNLM)
#'           - Intuition: MNLM as a set of binary logits (party id)
#'           - POE5 nels_small, Table 16.2 and 16.3
#'   Part 6  Conditional logit (CLM)
#'           - POE5 cola data, Table 16.4
#'           - Cracker data, Franses & Paap
#'           - Nested logit
#'           - Mixing individual- and alternative-specific variables
#'   Part 7  Ordered probit and logit (nels_small)
#'   Part 8  Censored data: the Tobit model (Example 16.16)
#'   Part 9  Extension: a Bayesian logit model (coke data)
#'
#' Run the script line by line. Each part starts with rm(list=ls()) and
#' reads its own data, so the parts can be run independently (but load
#' mosaic first, see below).
#' ===========================================================================

library(mosaic)   # tally(), favstats(), makeFun(), plotFun() and dplyr verbs

#' ===========================================================================
#' # Part 1: Binary choice -- the transport data
#' ===========================================================================
rm(list=ls())

#' Data definition file:
#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/transport.def")

# Obs:   21 
# 
# autotime	commute time via auto, minutes
# bustime	commute time via bus, minutes
# dtime		=(bus time - auto time)/10, 10 minute units
#			if positive, bus takes longer time than auto
#			if negative, bus is faster than auto
# auto		= 1 if auto chosen

#' Read the data
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/transport.rdata"))

#View(transport)
names(transport)
head(transport)

summary(transport)
tally(~auto, data=transport) # frequency
tally(~auto, data=transport, format="percent")
gf_histogram(~auto, data=transport)


#' ---------------------------------------------------------------------------
#' ## 1.1 Linear probability model (LPM)
#' ---------------------------------------------------------------------------
LPM <- lm(auto ~ dtime, data = transport)
summary(LPM)

round(confint(LPM),2)

#' We estimate that if travel times by public transportation and automobile
#' are equal, so that DTime =0, then the probability of a person choosing
#' automobile travel is 0.4848, close to 50-50, with 95% interval estimates
#' of [0.34, 0.63]
#' 
#' We estimate that, holding all else constant, an increase of 10 minutes
#' in the difference in travel time, increasing public transportation 
#' travel time relative to automobile travel time, increases the probability
#' of choosing automobile travel by 0.07, with a 95% interval estimates
#' of [0.04, 0.10], which seems relatively precise. 

#' The fitted model can be used to estimate the probability of 
#' automobile travel for any commuting time differentials. 
#' 
#' For Example, if dtime=0, no commuting time differentials 
f <- makeFun(LPM) 
f(dtime=0, interval="confidence")
f(0)

#' If dtime =1, a 10 minute longer commute by public transportation
#' we estimate the probability of automobile travel to be 0.5551. 
f(dtime=1)
f(1)

#' Exercise: what does the LPM predict for dtime = -10 and dtime = 10?
#' Why is this a problem?


#' ---------------------------------------------------------------------------
#' ## 1.2 Example 16.3: Probit maximum likelihood, a small example
#' ---------------------------------------------------------------------------
#' Before we use glm(), we estimate a probit model "by hand" on three
#' observations, to see what the software does behind the scenes.

# The data
y <- c(1,1,0)
x <- c(1.5,0.6,0.7)
one <- c(1,1,1)
X <- cbind(one,x)
X

# Log likelihood function requires a model matrix X and binary dependent variable y
probit.nll <- function(beta) {
  # linear predictor
  eta <- X %*% beta
  # probability
  p <- pnorm(eta)
  # negative log-likelihood
  -sum((1-y)*log(1-p)+y*log(p))
}

# The gradient requires a model matrix X and binary dependent variable y
probit.gr <- function(beta) {
  # linear predictor
  eta <- X %*% beta
  # probability
  p <- pnorm(eta)
  # chain rule
  u <- dnorm(eta)*(y-p)/(p*(1-p))
  # gradient
  -crossprod(X, u)
}

# glm estimate, remove default intercept
fit.probit <- glm(y ~ X - 1, family = binomial(link = "probit"))
summary(fit.probit)

# our own estimation using optim with parameter starting values b
b <- c(1,1)

fit <- optim(b, probit.nll, gr = probit.gr, method = "BFGS", hessian = TRUE)
fit

# comparison
unname(coef(fit.probit))
fit$par

# Variance cov matrix
vcov(fit.probit)
sqrt(diag(vcov(fit.probit)))

# The inverse of the Hessian of the negative log-likelihood is the
# estimated variance-covariance matrix
fisher_info <- solve(fit$hessian)
sqrt(diag(fisher_info))

rm(y, x, one, X, probit.nll, probit.gr, fit.probit, b, fit, fisher_info)


#' ---------------------------------------------------------------------------
#' ## 1.3 Probit model, transport data
#' ---------------------------------------------------------------------------
#?glm
p1 <- glm(auto ~ dtime, family = binomial(link = "probit"), data = transport)
summary(p1)

# x=TRUE stores the model matrix (X) in the object. The erer package below
# needs it, BUT marginaleffects returns zero slopes for a model fitted with
# x=TRUE. Therefore we keep two versions of the same model.
p1x <- glm(auto ~ dtime, x=TRUE, family = binomial(link = "probit"), data = transport)
head(p1x$x)

#' The negative sign of the intercept implies that when commuting times
#' by bus and auto are equal so that dtime =0, individuals have a bias 
#' against driving to work, relative to public transportation. 

#' The positive sign of b2 indicates that an increase in public 
#' transportation travel time, relative to auto travel time, increase 
#' the probability that an individual will choose to drive to work, and 
#' this coeff is statistically significant. 

confint(p1) # confidence interval on parameters

#' The estimated probability of a person choosing to drive to work when 
#' dtime = 0, POE5, p. 691: 
f <- makeFun(p1) 
f(dtime=0)
f(0)
f(-5:5)

# Predicted probability at the mean of dtime, compared with the sample share
f(mean(transport$dtime))
mean(transport$auto)

# Plot the data and the predicted probabilities
xyplot(auto~dtime, data=transport, scales=list(y=list(-1,1)))
plotFun(f(dtime)~dtime, add = TRUE)

#' Plot the predicted probabilities as a function of dtime
library(rockchalk)
predictOMatic(p1) # predicted probabilities
predictOMatic(p1, interval="confidence")

plotCurves(p1, plotx = "dtime")
plotCurves(p1, plotx = "dtime", opacity=80, col="red", interval="confidence") # Wald CI

# The same plot 
plotSlopes(p1, plotx = "dtime", interval = "conf")


#' ---------------------------------------------------------------------------
#' ## 1.4 Probit marginal effects
#' ---------------------------------------------------------------------------
#' Finding the marginal effect at dtime=2 by hand
#' (i.e., assuming travel via public transportation 
#' takes 20 minutes longer than auto travel)

f(2) 
qnorm(f(2)) # qnorm calculates the inverse of the cdf
dnorm(qnorm(f(2))) # then calculating the dnorm of the inverse gives you the pdf
dnorm(qnorm(f(2)))*coef(p1)[2] # finally, the marginal effect

#' A 10-minutes increase in the travel time via public
#' transportation increases the probability of travel via auto by
#' approximately 0.1037, given that taking the bus already requires 
#' 20 minutes more travel time than driving. 

# The marginal effect as a function of dtime
g <- function(x) {dnorm(qnorm(f(x)))*coef(p1)[2]}
g(2)
curve(g(x), -10,10, main="Plot of the marginal effect of dtime", xlab="dtime") 
segments(2,0,2,g(2), col="blue", lty=2)
segments(2,g(2),-10,g(2), col="red", lty=2)

#' ### The marginaleffects package
#browseURL("https://marginaleffects.com")
#browseURL("https://marginaleffects.com/chapters/slopes.html")

library(marginaleffects)
library(modelsummary)
modelsummary(p1) 

# The marginal effect at dtime=2, with standard error
slopes(p1, newdata = datagrid(dtime = 2))

# The slopes() function produces distinct estimates of the marginal effect
# for each row of the data used to fit the model.
#?slopes
mfx <- slopes(p1)
head(mfx)

# make a plot, add confidence intervals from conf.low and conf.high
mfx %>% 
  ggplot(aes(x = dtime, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2)

# The average marginal effect (AME) is the average of the marginal effects
# for each observation in the data set, POE5, p. 689
avg_slopes(p1)

tidy(mfx)
glance(mfx) # R2, AIC, BIC, etc.

# Marginal Effect at the Mean (MEM)
slopes(p1, newdata = "mean")
g(mean(transport$dtime)) # by hand

#' ### The same numbers with other packages
#' <http://www.rdocumentation.org/packages/erer/functions/maBina>
#' The function "maBina" calculates marginal effects for a binary 
#' probit or logit model and their standard errors.
#' x.mean =TRUE: calculate marginal effects at the means of independent variables.
#' If FALSE, marginal effects are calculated for each observation and then averaged. 
library(erer)
maBina(w = p1x, x.mean = TRUE, rev.dum = TRUE)  # calculated at the mean of a variable 
maBina(w = p1x, x.mean = FALSE, rev.dum = TRUE) # AME, POE5, p. 689 (POE4, p. 594), calculated as the mean of each data point

# erer::maTrend - Plot the predicted probabilities as a function of dtime
dp <- maBina(w = p1x, x.mean = TRUE, rev.dum = TRUE)
tdp <- maTrend(q = dp, nam.c = "dtime", simu.c = FALSE)
tdp
plot(tdp)

#' <http://www.rdocumentation.org/packages/mfx/functions/probitmfx>
#' The function "probitmfx" from the "mfx" package calculates 
#' the marginal effects for a probit regression.
library(mfx)
probitmfx(formula=auto ~ dtime, data=transport)  # calculated at the mean of a variable 
probitmfx(formula=auto ~ dtime, atmean = FALSE, data=transport) # AME, POE5, p. 689, calculated as the mean of each data point


#' ---------------------------------------------------------------------------
#' ## 1.5 Probit: marginal effects and elasticities with the delta method
#' ---------------------------------------------------------------------------
library(car)

# Set the specific value of dtime
dtime_value <- 2  # Replace with your desired value of dtime

# Use deltaMethod to calculate the marginal effect with standard error
# Define the parameter names as b1 (Intercept) and b2 (coefficient of dtime)
marginal_effect_probit <- deltaMethod(
  object = p1,
  g = paste("b2 * dnorm(b1 + b2 *", dtime_value, ")"),
  parameterNames = paste("b", 1:2, sep = "")
)

# The marginal effect with standard error
marginal_effect_probit

# marginaleffects package
slopes(p1, newdata = datagrid(dtime = 2))

# Use deltaMethod to calculate the elasticity with standard error
elasticity_probit <- 
  deltaMethod(
    object = p1,
    g = paste("(b2 * dnorm(b1 + b2 *", dtime_value, ")) *", dtime_value, "/ pnorm(b1 + b2 *", dtime_value, ")"),
    parameterNames = paste("b", 1:2, sep = ""))

# The elasticity with standard error
elasticity_probit

# marginaleffects package, slope = "eyex" gives the elasticity
slopes(p1, newdata = datagrid(dtime = 2), slope = "eyex")

# The probit elasticity as a function of dtime
ep <- function(x) {dnorm(qnorm(f(x)))*coef(p1)[2]*x/(f(x))}
plotFun(ep(dtime)~dtime, xlim = c(-10,10))
ep(-5:5)

#' Interpretation: the elasticity of 0.295 means that for a 1% increase in
#' dtime, the probability of using an automobile increases by 0.295%
#' relative to the current probability.
#' To calculate the new probability:
#' New Probability = Start Probability x (1 + Elasticity x Percentage Increase in dtime)

# Define the function to calculate the new probability
calculate_new_probability <- function(model, dtime_value, percentage_increase) {
  # Calculate the linear predictor at dtime_value
  lp_at_dtime <- coef(model)["(Intercept)"] + coef(model)["dtime"] * dtime_value
  # Calculate the current probability (CDF of the probit model)
  current_probability <- pnorm(lp_at_dtime)
  # Calculate the marginal effect at dtime_value (slope of the tangent)
  marginal_effect_at_dtime <- coef(model)["dtime"] * dnorm(lp_at_dtime)
  # Calculate the elasticity at dtime_value
  elasticity_at_dtime <- (marginal_effect_at_dtime * dtime_value) / current_probability
  # Calculate the new probability after the percentage increase in dtime
  new_probability <- current_probability * (1 + elasticity_at_dtime * (percentage_increase / 100))
  return(new_probability)
}

# Example usage
dtime_value <- 2           # dtime = 2 (20 minutes)
percentage_increase <- 10  # 10% increase in travel time

# Call the function for the probit model 'p1'
new_probability <- calculate_new_probability(p1, dtime_value, percentage_increase)

# Print the new probability
new_probability

# Starting point for the elasticity
f(dtime=2)

# increase in the probability for a 10% increase in travel time
new_probability-f(dtime=2)

#' ### Graphical view: the probit CDF and its tangent at dtime = 2
#' The slope of the tangent is the marginal effect.

# Define a sequence of dtime values for plotting the CDF
dtime_seq <- seq(-10, 10, length.out = 100)

# Compute the linear predictor (probit) at each value of dtime
lp <- coef(p1)["(Intercept)"] + coef(p1)["dtime"] * dtime_seq

# Compute the CDF (probit link) using pnorm
cdf <- pnorm(lp)

# Compute the linear predictor and the CDF at dtime_value = 2
lp_at_dtime <- coef(p1)["(Intercept)"] + coef(p1)["dtime"] * dtime_value
cdf_at_dtime <- pnorm(lp_at_dtime)

# Calculate the marginal effect at dtime_value = 2 (slope of the tangent)
marginal_effect_at_dtime <- coef(p1)["dtime"] * dnorm(lp_at_dtime)

# Create a data frame for ggplot
ggplot_df <- data.frame(dtime = dtime_seq, CDF = cdf)

# Generate the plot
ggplot(ggplot_df, aes(x = dtime, y = CDF)) +
  # Plot the CDF curve
  geom_line(linewidth = 1.2, color = "blue") +
  # Add a vertical line at dtime = 2
  geom_vline(xintercept = dtime_value, linetype = "dashed", color = "red") +
  # Add a point at the CDF value for dtime = 2
  annotate("point", x = dtime_value, y = cdf_at_dtime, size = 3, color = "red") +
  # Draw the tangent line at dtime = 2
  geom_abline(intercept = cdf_at_dtime - marginal_effect_at_dtime * dtime_value,
              slope = marginal_effect_at_dtime, linetype = "dotted", color = "black") +
  # Customize labels and title
  labs(x = "dtime", y = "Probit CDF",
       title = "CDF of Probit Model with Tangent at dtime = 2",
       subtitle = "Blue: CDF, Red: Vertical line, Black: Tangent line") +
  theme_minimal()


#' ---------------------------------------------------------------------------
#' ## 1.6 A closer look at the logistic function
#' ---------------------------------------------------------------------------
#' Some identities that are useful when working with logit elasticities.
#' With the logistic cdf L(z) = 1/(1+exp(-z)):  1 - L(z) = L(-z)

# logistic cdf, two equivalent ways of writing it
curve(exp(x)/(1+exp(x)),-2,2)
curve(1/(1+exp(-x)),-2,2, add=TRUE, col="red")

# logistic pdf, symmetric around zero
curve(exp(-x)/(1+exp(-x))^2,-2,2)
curve(exp(x)/(1+exp(x))^2,-2,2, add=TRUE, col="red")

# pdf = cdf(x)*(1-cdf(x))
curve((1/(1+exp(-x)))*(1-(1/(1+exp(-x)))),-2,2, add=TRUE, col="blue")

# Four ways of writing 1 - L(z)
p <- function(l) {1/(1+(1/exp(-l)))}
p(2)

g <- function(l) {exp(-l)/(1+exp(-l))}
g(2)

h <- function(l) {1/(1+exp(l))}
h(2)

i <- function(l) {1/(1+exp(-l))}
1-i(2)

x <- seq(-2,2, length.out = 100)
curve(p, -2,2)
lines(x,g(x), col="red")
lines(x,h(x), col="blue")
lines(x,1-i(x), col="green")

rm(x,p,g,h,i)  # remove x, g, h, i 

# Elasticity, with index z = b*x and b=1: pdf*b*x/cdf = b*x*(1-cdf)
curve((exp(x)/(1+exp(x))^2)*x*1/(1/(1+exp(-x))),-2,2)
curve(x*(1-(1/(1+exp(-x)))),-2,2, add=TRUE, col="red")

# Quasi elasticity (Franses & Paap): pdf*b*x = cdf*(1-cdf)*b*x
curve(x*(exp(x)/(1+exp(x))^2),-2,2, add=TRUE, col="blue")


#' ---------------------------------------------------------------------------
#' ## 1.7 Logit model, transport data
#' ---------------------------------------------------------------------------
l1 <- glm(auto ~ dtime, family = binomial(link = "logit"), data = transport)
summary(l1)

compareCoefs(p1,l1) # comparing the probit (model 1) and logit (model 2)
coef(l1)[2]/coef(p1)[2] # rule of thumb: logit coefficients are approx. 1.6 times probit

# Flip coding on auto (compare the estimates from logit models)
transport <- transport %>% mutate(inv.auto = ifelse(auto==0,1,0))
il1 <- glm(inv.auto ~ dtime, family = binomial(link = "logit"), data = transport)
summary(il1)
plotCurves(il1, plotx = "dtime", opacity=80, col="red", interval="confidence") # Wald CI

compareCoefs(l1,il1)

#' ### Logit marginal effects and elasticities
fl <- makeFun(l1) # predicted probability, logit

# Logit elasticity as a function of dtime: b*x*(1-P)
el <- function(x) { coef(l1)[2]*x*(1-fl(x)) }
plotFun(el(dtime)~dtime, xlim = c(-10,10))

dtime_value <- 2

# Use deltaMethod to calculate the marginal effect with standard error for logit
marginal_effect_logit <- deltaMethod(
  object = l1,
  g = paste("b2 * exp(b1 + b2 *", dtime_value, ") / (1 + exp(b1 + b2 *", dtime_value, "))^2"),
  parameterNames = paste("b", 1:2, sep = "")
)

# The marginal effect with standard error
marginal_effect_logit

# marginaleffects package
slopes(l1, newdata = datagrid(dtime = 2))

# Use deltaMethod to calculate the elasticity with standard error for logit
# Predicted probability P for the logit model
P_at_dtime_logit <- paste("1 / (1 + exp(-(b1 + b2 *", dtime_value, ")))")

# Define the formula for elasticity
elasticity_logit <- deltaMethod(
  object = l1,
  g = paste("(b2 * exp(b1 + b2 *", dtime_value, ") / (1 + exp(b1 + b2 *", dtime_value, "))^2) *", 
            dtime_value, "/ (", P_at_dtime_logit, ")"),
  parameterNames = paste("b", 1:2, sep = "")
)

# The elasticity with standard error
elasticity_logit
el(2) # by hand

# marginaleffects package
slopes(l1, newdata = datagrid(dtime = 2), slope = "eyex")

# Compare with the probit elasticity
elasticity_probit

# Average elasticity over all observations.
# Note that dtime takes both negative and positive values, so the
# individual elasticities change sign, and the average is hard to interpret.
avg_slopes(l1, slope = "eyex")

# The same by hand: elasticity at all observations, then averaged
mfx_l <- slopes(l1)
mean(mfx_l$estimate*transport$dtime/fl(dtime=transport$dtime))

# Note that the mean of dtime is negative, which makes the elasticity at this point negative
mean(transport$dtime)

# The marginal effect at the mean
dydx <- slopes(l1, newdata = "mean")
dydx$estimate

# Predicted probability at the mean
fl(dtime=mean(transport$dtime))

# elasticity at mean
dydx$estimate*mean(transport$dtime)/fl(dtime=mean(transport$dtime))


#' ===========================================================================
#' # Part 2: Example 16.6, choice between Coke and Pepsi
#' ===========================================================================
rm(list=ls())

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/coke.def")

# Obs:   1140 individuals
# 
# coke       	=1 if coke chosen, =0 if pepsi chosen
# pr_pepsi        price of 2 liter bottle of pepsi
# pr_coke         price of 2 liter bottle of coke
# disp_pepsi      = 1 if pepsi is displayed at time of purchase, otherwise = 0
# disp_coke       = 1 if coke is displayed at time of purchase, otherwise = 0
# pratio          price coke relative to price pepsi

#' Read the data
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/coke.rdata"))
#View(coke)

#logit <- glm(coke~I(pr_coke/pr_pepsi)+disp_coke+disp_pepsi, family = binomial(link = "logit"), data = coke)
logit <- glm(coke~pratio+disp_coke+disp_pepsi, family = binomial(link = "logit"), data = coke)
summary(logit)

# The same model with x=TRUE, needed by erer::maBina below
# (do not use this version with marginaleffects, it gives zero slopes)
logitx <- glm(coke~pratio+disp_coke+disp_pepsi, family = binomial(link = "logit"), data = coke, x=TRUE)

#' Here we do not interpret the coefficients, or the magnitude of the coeff. 
#' We can only interpret the sign of the coeff., i.e. whether the explanatory 
#' variables make the choice more likely or less likely.

# CIs using profiled log-likelihood
confint(logit)
# CIs using standard errors
confint.default(logit)


#' ---------------------------------------------------------------------------
#' ## 2.1 Hypothesis tests
#' ---------------------------------------------------------------------------
# Wald hypothesis test, POE5, p. 695
# <http://www.rdocumentation.org/packages/aod/functions/wald.test>
# Computes a Wald chi-square test for 1 or more coefficients, given their variance-covariance matrix.
library(aod)
wald.test(b = coef(logit), Sigma = vcov(logit), Terms = 2:4) # overall significance (except intercept)

#find Chi-Square critical value
qchisq(0.95, df=3, lower.tail=TRUE)

#' Since x2 > CV, we reject the null hypothesis that none of the 
#' explanatory variables help explain the choice of Coke versus Pepsi. 

probit <- glm(coke~pratio+disp_coke+disp_pepsi, family = binomial(link = "probit"), data = coke)
summary(probit)

#' Example 16.7: Test the hypothesis disp_coke + disp_pepsi = 0, i.e. that the
#' Coke and Pepsi displays have equal but opposite effects, POE5, p. 696
library(multcomp)
summary(glht(probit, linfct = ("disp_coke+disp_pepsi = 0")))

#' Joint hypothesis, disp_coke=0 & disp_pepsi=0
#' Test the joint hypothesis that neither the Coke nor Pepsi
#' display affects the probability of choosing Coke.
wald.test(b = coef(probit), Sigma = vcov(probit), Terms = 3:4) 
#' We reject H0 and conclude that the Coke and Pepsi displays have an effect on 
#' the probability of choosing Coke. 


#' ---------------------------------------------------------------------------
#' ## 2.2 Results from both models
#' ---------------------------------------------------------------------------
library(stargazer)
stargazer(probit,logit, type = "text", intercept.bottom = FALSE)
#' The parameters and their estimates vary across the models and no 
#' direct comparison is very useful, but some rules of thumb exist,
#' See the PPt (page = 29). 


#' ---------------------------------------------------------------------------
#' ## 2.3 Predicted probabilities
#' ---------------------------------------------------------------------------
library(rockchalk)
mean(coke$pratio)
hist(coke$pratio)

plotCurves(logit, plotx = "pratio")
plotCurves(logit, plotx = "pratio", opacity=80, col="red", interval="confidence") # Wald CI

plotCurves(logit, plotx = "pratio", modx = "disp_pepsi")
abline(v=mean(coke$pratio), col="red") # the mean of pratio

plotCurves(logit, plotx = "pratio", modx = "disp_coke")
abline(v=mean(coke$pratio), col="red")

# Predicted probabilities, effect of a Pepsi display at the mean of pratio
f <- makeFun(logit) 
f(pratio=mean(coke$pratio), disp_pepsi = 1, disp_coke = mean(coke$disp_coke))
f(pratio=mean(coke$pratio), disp_pepsi = 0, disp_coke = mean(coke$disp_coke))

P1=f(pratio=mean(coke$pratio), disp_pepsi = 1, disp_coke = mean(coke$disp_coke))
P0=f(pratio=mean(coke$pratio), disp_pepsi = 0, disp_coke = mean(coke$disp_coke))

P1-P0  

plotCurves(logit, plotx = "pratio", modx = "disp_pepsi")
abline(v=mean(coke$pratio))
abline(h=P1, col="red")
abline(h=P0, col="blue")

# Effect of a discount: pratio=0.8 vs pratio=0.9, with a Coke display
P1=f(pratio=0.8, disp_pepsi = 0, disp_coke = 1)
P0=f(pratio=0.9, disp_pepsi = 0, disp_coke = 1)
P1-P0

P1=f(pratio=0.8, disp_pepsi = mean(coke$disp_pepsi), disp_coke = 1)
P0=f(pratio=0.9, disp_pepsi = mean(coke$disp_pepsi), disp_coke = 1)
P1-P0

#' Exercise: What is the effect of disp_pepsi when pratio is 1.5?


#' ---------------------------------------------------------------------------
#' ## 2.4 Marginal effects
#' ---------------------------------------------------------------------------
#' <http://www.rdocumentation.org/packages/erer/functions/maBina>
#' x.mean =TRUE: calculate marginal effects at the means of independent variables.
#' If FALSE, marginal effects are calculated for each observation and then averaged.
library(erer)
maBina(w = logitx, x.mean = TRUE, rev.dum = TRUE)

#' The functions "probitmfx" and "logitmfx" from the "mfx" package calculate 
#' the marginal effects for probit and logit regressions.
#' <http://www.rdocumentation.org/packages/mfx/functions/probitmfx>
library(mfx)
logitmfx(formula=coke~pratio+disp_pepsi+disp_coke, data=coke)

library(marginaleffects)
avg_slopes(logit)                 # AME
slopes(logit, newdata = "mean")   # MEM


#' ---------------------------------------------------------------------------
#' ## 2.5 Elasticities in the logit model
#' ---------------------------------------------------------------------------
#' For a logit model, the elasticity of P with respect to x is b*x*(1-P)
logit

#' Elasticity as a function of pratio
#' Note: 1-P = 1/(1+exp(index)), see section 1.6
e <- function(pratio,disp_coke,disp_pepsi) {
  1/(1+exp(coef(logit)[1]+pratio*coef(logit)[2]+disp_coke*coef(logit)[3]+
             disp_pepsi*coef(logit)[4]))*pratio*coef(logit)[2]}
summary(coke$pratio)
#' pratio=(pr_coke/pr_pepsi)

#' No advertising
curve(e(x,0,0), min(coke$pratio), max(coke$pratio), ylim=c(-4,0), xlab = "pratio", ylab = "elasticity",
      main = "Plot of the pratio elasticity")

#' Coke advertising
curve(e(x,1,0), min(coke$pratio), max(coke$pratio), add=T, col="blue") # Own advertising make choice of coke more inelastic

#' Pepsi advertising
curve(e(x,0,1), min(coke$pratio), max(coke$pratio), add=T, col="red") # Competitor advertising make choice of coke more elastic

#' Both advertise
curve(e(x,1,1), min(coke$pratio), max(coke$pratio), add=T, col="green")
legend("topright", legend=c("No advertising", "Coke advertising", "Pepsi advertising","Both advertise"),
       col=c("black","blue","red", "green"), lty=c(1,1,1,1))

#' At equal price, pepsi advertising makes choice of coke more elastic
e(1,1,0) #pratio =1 (i.e., pr_coke = pr_pepsi)
e(1,1,1)

#' Elasticity as a function of pratio, ver. 2, using the predicted probability
e2 <- function(pratio,disp_coke,disp_pepsi) {(1-f(pratio,disp_coke,disp_pepsi))*pratio*coef(logit)[2]}
e2(1,1,0)
e2(1,1,1)


#' ---------------------------------------------------------------------------
#' ## 2.6 A logit model with individual prices: own and cross price elasticities
#' ---------------------------------------------------------------------------
names(coke)
logit2 <- glm(coke ~ pr_coke + pr_pepsi + disp_coke + disp_pepsi, family = binomial(link = "logit"), data = coke)
summary(logit2)

g <- makeFun(logit2)

#' Own price elasticity
e3 <- function(pr_coke,pr_pepsi,disp_coke,disp_pepsi) {(1-g(pr_coke,pr_pepsi,disp_coke,disp_pepsi))*pr_coke*coef(logit2)[2]}

curve(e3(x,mean(coke$pr_pepsi),0,0), min(coke$pr_coke), max(coke$pr_coke), ylim=c(-2,0), xlab = "pr_coke", main="Coke own price elasticity") 
curve(e3(x,mean(coke$pr_pepsi),1,0), min(coke$pr_coke), max(coke$pr_coke), add=T, col="red") # Own advertising make choice of coke more inelastic

#' Substitute (cross) price elasticity
e4 <- function(pr_coke,pr_pepsi,disp_coke,disp_pepsi) {(1-g(pr_coke,pr_pepsi,disp_coke,disp_pepsi))*pr_pepsi*coef(logit2)[3]}
curve(e4(mean(coke$pr_coke),x,0,0), min(coke$pr_pepsi), max(coke$pr_pepsi), ylim=c(0,2), xlab = "pr_pepsi", main="Coke: substitute pepsi price elasticity") 

#' ### The same with the marginaleffects package
# AME
avg_slopes(logit2)
plot_predictions(logit2, condition = "pr_coke")

# Average elasticities
avg_slopes(logit2, slope = "eyex")

#https://marginaleffects.com/articles/slopes.html#marginal-effect-at-user-specified-values
# Here only pr_coke varies, the values of the other variables are shown in the output
slopes(logit2,
       slope = "eyex",
       newdata = datagrid(pr_coke = seq(0.8,1.8, by=0.1))) # expand with the other variables

# Elasticities, plotted against the own price
plot_slopes(logit2, slope = "eyex", variables = "pr_coke", condition = "pr_coke")

plot_slopes(logit2, slope = "eyex", variables = "pr_pepsi", condition = "pr_pepsi")

# Own price elasticity of coke, as a function of the price of pepsi,
# with and without displays.
# This is a ggplot object, hence we can fix axes and labels
plot_slopes(logit2, slope = "eyex", variables = "pr_coke",
            condition = list("pr_pepsi" = seq(0.8,1.8, by=0.1), "disp_coke" = 0:1, "disp_pepsi" = 0:1)) +
  ylim(-1.6,0) + ylab("Elasticity") + xlab("Price of Pepsi")

#' Exercise: Make a plot of the own price elasticity of pepsi, from the choice
#' of pepsi as the dependent variable.
#' What is the own price elasticity at the mean price of pepsi, with no advertising?


#' ===========================================================================
#' # Part 3: Binary choice in marketing, Heinz vs Hunts ketchup
#' ===========================================================================
#' Franses & Paap, Chapter 4 (section 4.4)
rm(list=ls())

library(broom)
library(tidyverse)

url <- "https://raw.githubusercontent.com/oysteinm/data/refs/heads/main/heinzhunts.csv"
heinzhuntsdata <- read_csv(url, id=NULL)

head(heinzhuntsdata)
names(heinzhuntsdata)
table(heinzhuntsdata$HOUSEHOLDID)
histogram(~HOUSEHOLDID, data=heinzhuntsdata, width=1, type="count")

# The last purchase of each household is kept as a hold-out sample
tally(~LASTPURCHASE, data=heinzhuntsdata)
df.within <- heinzhuntsdata %>% filter(LASTPURCHASE == 0)

#' Logit model
mod.logit <- glm(HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+
                   FEATDISPLHUNTS+log(PRICEHEINZ/PRICEHUNTS),
                 x=TRUE, family = binomial(link = "logit"), data = df.within)

summary(mod.logit)
glance(mod.logit)

# CIs using profiled log-likelihood
confint(mod.logit) # confidence interval on parameters

# CIs using standard errors
confint.default(mod.logit)

#' Probit model
mod.probit <- glm(HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+
                    FEATHUNTS+FEATDISPLHUNTS+log(PRICEHEINZ/PRICEHUNTS),
                  x=TRUE, family = binomial(link = "probit"), data = df.within)

summary(mod.probit)
glance(mod.probit)

confint(mod.probit)

# Are prices lower when Heinz is on display?
summary(lm(PRICEHEINZ~DISPLHEINZ, data = df.within))

#' Marginal effects
#' <http://www.rdocumentation.org/packages/erer/functions/maBina>
library(erer)
maBina(w = mod.logit, x.mean = TRUE, rev.dum = TRUE)  # calculated at the mean of a variable 
maBina(w = mod.logit, x.mean = FALSE, rev.dum = TRUE) # AME, calculated as the mean of each data point

histogram(~log(PRICEHEINZ/PRICEHUNTS), width=0.01, data=heinzhuntsdata, type="count")

maBina(w = mod.probit, x.mean = TRUE, rev.dum = TRUE)  # calculated at the mean of a variable 
maBina(w = mod.probit, x.mean = FALSE, rev.dum = TRUE) # AME, calculated as the mean of each data point

#' <http://www.rdocumentation.org/packages/mfx/functions/probitmfx>
library(mfx)
probitmfx(formula=HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+log(PRICEHEINZ/PRICEHUNTS),
          data=df.within)  # calculated at the mean of a variable 
probitmfx(formula=HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+log(PRICEHEINZ/PRICEHUNTS),
          atmean = FALSE, data=df.within) # AME, calculated as the mean of each data point

#' To plot against the log price ratio, we make it a variable in the data
df.within %>% mutate(lnP_HEINZ_HUNTS=log(PRICEHEINZ/PRICEHUNTS)) -> df.within

fit <- glm(HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+lnP_HEINZ_HUNTS,
           x=TRUE, family = binomial(link = "logit"), data = df.within)
summary(fit)

#' Plot the predicted probabilities as a function of the log price ratio
library(rockchalk)
predictOMatic(fit) # predicted probabilities
predictOMatic(fit, interval="confidence")

plotCurves(fit, plotx = "lnP_HEINZ_HUNTS")
plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", opacity=80, col="red", interval="confidence") # Wald CI

plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="DISPLHEINZ", xlim=c(-0.5,1.5))
plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="FEATHEINZ", xlim=c(-0.5,1.5))
plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="FEATDISPLHEINZ", xlim=c(-0.5,1.5))

plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="DISPLHUNTS", xlim=c(-0.5,1.5))
plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="FEATHUNTS", xlim=c(-0.5,1.5))
plotCurves(fit, plotx = "lnP_HEINZ_HUNTS", modx ="FEATDISPLHUNTS", xlim=c(-0.5,1.5))

plotModel(fit, HEINZ~lnP_HEINZ_HUNTS, system="g")

#' Predicted probabilities
f <- makeFun(fit) # mosaic
f(DISPLHEINZ=1, FEATHEINZ=0, FEATDISPLHEINZ=0, DISPLHUNTS=0, FEATHUNTS=0, FEATDISPLHUNTS=0, lnP_HEINZ_HUNTS=0)
f(DISPLHEINZ=0, FEATHEINZ=0, FEATDISPLHEINZ=0, DISPLHUNTS=0, FEATHUNTS=0, FEATDISPLHUNTS=0, lnP_HEINZ_HUNTS=0)
f(0,0,0,0,0,0,0)

summary(df.within$lnP_HEINZ_HUNTS)

#' Elasticities
#' NB: the regressor is already in logs, x = ln(PRICEHEINZ/PRICEHUNTS).
#' Then dP/dx = b*P*(1-P), and the elasticity of P with respect to the
#' price RATIO is simply b*(1-P). The curves below multiply by x as well,
#' i.e. they are the elasticity with respect to the LOG price ratio.
#' Discuss: which one is the economically meaningful elasticity?
curve((1-f(0,0,0,0,0,0,x))*x*coef(fit)[8], min(df.within$lnP_HEINZ_HUNTS), max(df.within$lnP_HEINZ_HUNTS)) # plot of the elasticity

e <- function(lnP_HEINZ_HUNTS)
{(1-f(0,0,0,0,0,0,lnP_HEINZ_HUNTS))*lnP_HEINZ_HUNTS*coef(fit)[8]}

curve(e(x), min(df.within$lnP_HEINZ_HUNTS), max(df.within$lnP_HEINZ_HUNTS)) # plot of the elasticity

# The same, writing out 1-P = 1/(1+exp(index))
g <- function(DISPLHEINZ, FEATHEINZ, FEATDISPLHEINZ, DISPLHUNTS, FEATHUNTS, FEATDISPLHUNTS, lnP_HEINZ_HUNTS) 
{1/(1+exp(coef(fit)[1]+DISPLHEINZ*coef(fit)[2]+FEATHEINZ*coef(fit)[3]+FEATDISPLHEINZ*coef(fit)[4]+
            DISPLHUNTS*coef(fit)[5]+FEATHUNTS*coef(fit)[6]+FEATDISPLHUNTS*coef(fit)[7]+lnP_HEINZ_HUNTS*coef(fit)[8]))*lnP_HEINZ_HUNTS*coef(fit)[8]}

curve(g(0,0,0,0,0,0,x), min(df.within$lnP_HEINZ_HUNTS), max(df.within$lnP_HEINZ_HUNTS)) # plot of the elasticity

# h = P(y=0)
h <- function(DISPLHEINZ, FEATHEINZ, FEATDISPLHEINZ, DISPLHUNTS, FEATHUNTS, FEATDISPLHUNTS, lnP_HEINZ_HUNTS) 
{1/(1+exp(coef(fit)[1]+DISPLHEINZ*coef(fit)[2]+FEATHEINZ*coef(fit)[3]+FEATDISPLHEINZ*coef(fit)[4]+
            DISPLHUNTS*coef(fit)[5]+FEATHUNTS*coef(fit)[6]+FEATDISPLHUNTS*coef(fit)[7]+lnP_HEINZ_HUNTS*coef(fit)[8]))}

f(0,0,0,0,0,0,2)   # P(y=1)
1-h(0,0,0,0,0,0,2) # P(y=1) = 1 - P(y=0)

# Franses & Paap, p. 70: Quasi elasticity, P*(1-P)*b*x  (see the NB above)
curve(f(0,0,0,0,0,0,x)*(1-f(0,0,0,0,0,0,x))*x*coef(fit)[8], min(df.within$lnP_HEINZ_HUNTS), max(df.within$lnP_HEINZ_HUNTS), xlim=c(-0.5,1.5)) # plot of the elasticity


#' ---------------------------------------------------------------------------
#' ## 3.1 Model evaluation
#' ---------------------------------------------------------------------------
#browseURL("https://www.r-bloggers.com/evaluating-logistic-regression-models/")

# Pseudo R2 measures
library(pscl)
pR2(fit)

# Confusion matrix, with cutoff 0.5 and with cutoff 0.891,
# approximately the sample share of Heinz purchases
mean(df.within$HEINZ)
table(df.within$HEINZ, predict(fit, type = 'response') > 0.5)
table(df.within$HEINZ, predict(fit, type = 'response') > 0.891)

# Accuracy within sample
accuracy <- table(predict(fit, type = 'response') > 0.891, df.within$HEINZ)
sum(diag(accuracy))/sum(accuracy)

# Compare w cutoff =0.5
accuracy.default <- table(predict(fit, type = 'response') > 0.5, df.within$HEINZ)
sum(diag(accuracy.default))/sum(accuracy.default)

# Accuracy in the hold-out sample (the last purchase of each household)
df.holdout <- heinzhuntsdata %>% filter(LASTPURCHASE == 1)
df.holdout %>% mutate(lnP_HEINZ_HUNTS=log(PRICEHEINZ/PRICEHUNTS)) -> df.holdout

accuracy.holdout <- table(predict(fit, type = 'response', newdata = df.holdout) > 0.891, df.holdout$HEINZ)
sum(diag(accuracy.holdout))/sum(accuracy.holdout)

# ROC Curve
library(ROCR)
ROCRpred <- prediction(predict(fit, type = 'response') , df.within$HEINZ)
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))

#' What happens if we change the dep var to HUNTS?
tidy(glm(HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+lnP_HEINZ_HUNTS,
         family = binomial(link = "logit"), data = df.within))
tidy(glm(HUNTS~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+lnP_HEINZ_HUNTS,
         family = binomial(link = "logit"), data = df.within))

glance(glm(HEINZ~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+lnP_HEINZ_HUNTS,
           family = binomial(link = "logit"), data = df.within))
glance(glm(HUNTS~DISPLHEINZ+FEATHEINZ+FEATDISPLHEINZ+DISPLHUNTS+FEATHUNTS+FEATDISPLHUNTS+lnP_HEINZ_HUNTS,
           family = binomial(link = "logit"), data = df.within))


#' ===========================================================================
#' # Part 4: Interaction terms in a logit model
#' ===========================================================================
#' The Cracker data from the mlogit package. We keep the purchases of
#' Sunshine and Keebler only, and model the choice of Sunshine.
rm(list=ls())

library(mlogit)
data("Cracker", package = "mlogit")

head(Cracker)

with(Cracker, table(choice))

cr <- Cracker %>%
  filter(choice %in% c("sunshine", "keebler")) %>%
  mutate(choice2 = ifelse(choice=="sunshine",1,0))

logitMod <- glm(choice2 ~ price.sunshine*disp.sunshine, data=cr, family=binomial(link="logit"))
summary(logitMod)

#' In a nonlinear model the coefficient on the interaction term is not the
#' "interaction effect" on the probability. Look at the marginal effect of
#' price, with and without a display:
library(marginaleffects)
avg_slopes(logitMod, variables = "price.sunshine", by = "disp.sunshine")
plot_predictions(logitMod, condition = c("price.sunshine", "disp.sunshine"))


#' ===========================================================================
#' # Part 5: Multinomial logit models - MNLM
#' ===========================================================================

#' ---------------------------------------------------------------------------
#' ## 5.1 Intuition: the MNLM as a set of binary logits (party identification)
#' ---------------------------------------------------------------------------
rm(list=ls())

# Install haven if you haven't already:
# install.packages("haven")
library(haven)

# Use the raw URL
url <- "https://raw.githubusercontent.com/oysteinm/data/main/partyid4.dta"

# Read directly into a data frame / tibble
partyid4 <- read_dta(url)

# View the data
head(partyid4)
names(partyid4)

# All 3 groups
str(partyid4$party3)
tally(~party3, data = partyid4)
tally(~party3, data = partyid4, format = "percent")

# Democrats | Independent
tally(~dem_ind, data = partyid4)
# Excludes 538 Republicans
summary(glm(dem_ind ~ income, data = partyid4, family = binomial(link = "logit")))

# Republicans | Independent
tally(~rep_ind, data = partyid4)
# Excludes 693 Democrats
summary(glm(rep_ind ~ income, data = partyid4, family = binomial(link = "logit")))

# Democrats | Republicans
tally(~dem_rep, data = partyid4)
# Excludes 151 Independents
summary(glm(dem_rep ~ income, data = partyid4, family = binomial(link = "logit")))

# Fitting a MNLM
library(nnet)
summary(multinom(party3 ~ income, data=partyid4)) # uses the original data

# Alternative
library(mlogit)

# For data in wide format, mlogit names the alternatives after the values
# of the choice variable. So we first give party3 labels, in a copy of the data.
# Check that the counts match the tally above: 693 Dem, 151 Ind, 538 Rep.
partyid4_ml <- partyid4 %>% 
  mutate(party3 = factor(party3, levels = 1:3, labels = c("Democrat","Independent","Republican")))
tally(~party3, data = partyid4_ml)

DF <- mlogit.data(partyid4_ml, choice = "party3", shape="wide")
levels(DF$alt)

# Note the data structure, long or stacked, with 3 choices/categories per respondent
head(DF)

# The formula: no alternative-specific variables (0), income is individual specific
# Can choose any reference level
mod_dem <- mlogit(party3 ~ 0 | income, data = DF, reflevel = "Democrat")
summary(mod_dem)
head(model.matrix(mod_dem)) # the design matrix, J-1 intercepts and income slopes
summary(mlogit(party3 ~ 0 | income, data = DF, reflevel = "Republican"))
summary(mlogit(party3 ~ 0 | income, data = DF, reflevel = "Independent"))

# Reflevel Independent, compare MNLM with logit
# Democrat:income        -0.0027241  0.0037162 -0.7330    0.4635    
summary(glm(dem_ind ~ income, data = partyid4, family = binomial(link = "logit")))


#' ---------------------------------------------------------------------------
#' ## 5.2 POE5: choice of postsecondary education, nels_small
#' ---------------------------------------------------------------------------
rm(list=ls())

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/nels_small.def")

# Obs:   1000 observations 
# 
# psechoice	= 1 if first postsecondary education was no college
#           = 2 if first postsecondary education was a 2-year college
#           = 3 if first postsecondary education was a 4-year college
# hscath		= 1 if catholic high school graduate
# grades		= average grade in math, english and social studies on 13 point scale with 1 = highest
# faminc		= gross 1991 family income (in $1000)
# famsiz		= number of family members
# parcoll		= 1 if most educated parent graduated from college or had an advanced degree
# female		= 1 if female
# black		  = 1 if black

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/nels_small.rdata"))

# Note the data structure, 1 obs per respondent
head(nels_small)

#' Descriptive statistics
tally(~psechoice, data=nels_small)
tally(~psechoice, data=nels_small, format = "percent")
#' Of 1000 respondents, 22.2% selected not to attend a college upon graduation,
#'  25.1% selected to attend a 2-year college and 52.7% selected a 4-year college. 

summary(nels_small)
mean(nels_small$grades)
sd(nels_small$grades)

#' Two different ways to estimate the multinomial logit model 

#' ### First approach: using the mlogit package 
library(mlogit)

# For data in wide format, mlogit names the alternatives after the values
# of the choice variable. So we give psechoice labels, in a copy of the data
# (we keep psechoice numeric in nels_small for later use).
nels_ml <- nels_small %>% 
  mutate(psechoice = factor(psechoice, levels = 1:3,
                            labels = c("no college","two-year college","four-year college")))

DF <- mlogit.data(nels_ml, choice = "psechoice", shape="wide")
levels(DF$alt)

# Note the data structure, long or stacked, with 3 choices per respondent
head(DF,12)

# grades are individual specific
# grades = 1 A+
# grades = 13 F
f1 <- psechoice ~ 0 | grades

m1 <- mlogit(f1, data = DF, reflevel = "no college")
summary(m1)  # Table 16.2, page 706
head(model.matrix(m1)) # the design matrix

#' Based on the estimates, what can we say?
#' Recall that a larger numerical value of GRADES represents a poorer 
#' academic performance. 
#' The parameter estimates for the coefficients of GRADES are negative and 
#' statistically significant. 
#' 
#' If the value of GRADES increases, the probability that high-school
#' graduates will choose a 2-year or a 4-year college goes down, relative
#' to the probability of not attending college. 
#' This is the anticipated effect, as we expect that a poorer academic 
#' performance will increase the odds of not attending college. 

#' ### The choice of reference category
m2 <- mlogit(f1,DF, reflevel = "two-year college")
summary(m2)

m3 <- mlogit(f1,DF, reflevel = "four-year college")
summary(m3)

# The coefficients are log odds of an alternative relative to the
# reference category, ln(Pj/Pref), i.e. the reference is in the denominator
coef(m1) # ref. no college
coef(m1)[1] # two year/no
coef(m1)[2] # four year/no

coef(m2) # ref. two-year
coef(m2)[1] # no/two year
coef(m2)[2] # four year/two year

coef(m3) # ref. four year
coef(m3)[1] # no/four year
coef(m3)[2] # two year/four year

# The estimates are consistent across reference categories:
# ln(P2/P1) + ln(P3/P2) = ln(P3/P1)
# intercepts
all.equal(as.numeric(coef(m1)[1] + coef(m2)[2]), as.numeric(coef(m1)[2]))
# slopes
coef(m1)[3] + coef(m2)[4]
coef(m1)[4]
all.equal(as.numeric(coef(m1)[3] + coef(m2)[4]), as.numeric(coef(m1)[4]))
# etc.

#' ### Predicted probabilities
tally(~psechoice, format = "percent", data=nels_small) # actual
apply(fitted(m1, outcome=FALSE),2,mean) # predicted, average over individuals

# predicted values for each individual
head(fitted(m1))                 # probability of the chosen alternative
head(fitted(m1, outcome=FALSE))  # probabilities of all alternatives

#' ### Marginal effects, formula 16.20
#' The data for prediction has one row per alternative. Since grades is
#' individual specific, it takes the same value in all three rows.

# At the mean grade
z <- data.frame(alt=1:3, grades=rep(mean(nels_small$grades),3))
z
effects(m1, covariate="grades", data=z)
sum(effects(m1, covariate="grades", data=z)) # the marginal effects sum to zero
predict(m1, data=z) # probabilities at z (mean grades)

# The same marginal effects by hand, formula 16.20
predict(m1, data=z)[1]*(0-(0*predict(m1, data=z)[1]+coef(m1)[3]*predict(m1, data=z)[2]+
                             coef(m1)[4]*predict(m1, data=z)[3]))
predict(m1, data=z)[2]*(coef(m1)[3]-(0*predict(m1, data=z)[1]+coef(m1)[3]*predict(m1, data=z)[2]+
                                       coef(m1)[4]*predict(m1, data=z)[3]))
predict(m1, data=z)[3]*(coef(m1)[4]-(0*predict(m1, data=z)[1]+coef(m1)[3]*predict(m1, data=z)[2]+
                                       coef(m1)[4]*predict(m1, data=z)[3]))

# At the median grade (Table 16.3), grades = 6.64
z1 <- data.frame(alt=1:3, grades=rep(6.64,3))
effects(m1, covariate="grades", data=z1) 
predict(m1, data=z1)

# At the 5th percentile (Table 16.3), grades = 2.635
z2 <- data.frame(alt=1:3, grades=rep(2.635,3))
effects(m1, covariate="grades", data=z2) 
predict(m1, data=z2)

# Compare with the sample percentiles, slightly different due to how they are calculated
median(nels_small$grades)
quantile(nels_small$grades, c(0.05, 0.95))

#' ### Average marginal effects (AME) with standard errors
# https://stackoverflow.com/questions/54079553/how-to-get-average-marginal-effects-ames-with-standard-errors-of-a-multinomial

# marginal effects at all data points, then averaged
head(effects(m1, covariate="grades", data=DF))
colMeans(effects(m1, covariate="grades", data=DF))

# vector of covariate names
c.names <- "grades"

# get marginal effects for each id
ME.mnl <- sapply(c.names, function(x) 
  stats::effects(m1, covariate=x, data=DF), 
  simplify=FALSE) 

# get AMEs
(AME.mnl <- t(sapply(ME.mnl, colMeans)))

# standard errors by the delta method, using a numerical derivative
AME.fun <- function(betas) {
  tmp <- m1
  dframe <- DF
  tmp$coefficients <- betas
  ME.mnl <- sapply(c.names, function(x) 
    effects(tmp, covariate = x, data = dframe), simplify = FALSE)
  c(sapply(ME.mnl, colMeans))
}

library(numDeriv)
grad <- jacobian(AME.fun, m1$coef)
(AME.mnl.se <- matrix(sqrt(diag(grad %*% vcov(m1) %*% t(grad))), nrow = 1, byrow = TRUE))
AME.mnl / AME.mnl.se # z-values


#' ### Second approach: using the "multinom" function from the "nnet" package 
library(nnet)

nels_small$psechoice <- relevel(as.factor(nels_small$psechoice), ref = "1") #base/reference category
str(nels_small)

mod1 <- multinom(psechoice ~ grades, data = nels_small) # uses the original data!
summary(mod1)

# change basis category to 2-year college, compare with mlogit m2 above
nels_small$psechoice2 <- relevel(factor(nels_small$psechoice), ref = "2")
mod2 <- multinom(psechoice2~grades, data=nels_small) 
summary(mod2)
summary(m2)

#' A nice presentation 
library(stargazer)
stargazer(mod1, type = "text", intercept.bottom = FALSE,title = "The Multinomial nels_small Choice Estimates", header=FALSE)

# The multinom function does not include p-values for the regression coefficients,
# we can calculate p-values using Wald tests (here z-tests)
z <- summary(mod1)$coefficients/summary(mod1)$standard.errors
z

# 2-tailed p-values
p <- 2*(1 - pnorm(abs(z)))
p

#' MNLM with only the intercept included in the model
mod0 <- multinom(psechoice ~ 1, data = nels_small)
summary(mod0)

# Overall significance of the model, likelihood ratio test
anova(mod1,mod0)

#' Prediction of probabilities
# Use the function "fitted" to predict the probabilities within the data.
head(fitted(mod1)) # predicted probabilities per individual, slightly different from mlogit

predict(mod1, newdata=data.frame(grades=6.64), "probs") 
predict(mod1, newdata=data.frame(grades=2.635), "probs") # same as in Table 16.3

# At the same time 
predict(mod1, newdata=data.frame(grades=c(6.64,2.635)), type="probs")

# prediction at the mean
predict(mod1, newdata=data.frame(grades=mean(nels_small$grades)), type="probs")

#' Predicted probabilities as a function of grades
#browseURL("http://www.jstor.org/stable/pdf/25046697.pdf")
#browseURL("https://cran.r-project.org/web/packages/effects/")
library(effects)
plot(effect("grades", mod1, xlevels=list(grades=1:13)), lty=4, colors = "red",
     main = "Effect of a change in grades on the choice probabilities")

summary(Effect("grades", mod1, confint=TRUE))

# Same as Table 16.3
eff.models2 <- allEffects(mod1, xlevels=list(grades=c(6.64, 2.635)))
eff.models2
as.data.frame(eff.models2[[1]])

#' The probability ratio, or relative risk ratio
exp(coef(mod1))
#' The relative risk ratio for a one-unit increase in the variable 
#' grades is 0.7344448 for being in two-year college vs. no college.
#' The relative risk ratio for a one-unit increase in the variable grades 
#' is 0.4934929 for being in four-year college vs. no college.

#' Estimate many categories fast:
#browseURL("https://cran.r-project.org/web/packages/mnlogit/vignettes/mnlogit.pdf")


#' ===========================================================================
#' # Part 6: Conditional logit - CLM
#' ===========================================================================
#' In the MNLM, each explanatory variable has a different effect on each outcome.
#' And each explanatory variable (x) was individual specific.
#' Hence, in the MNLM there are J-1 parameters (b) for each x, but only a
#' single value of x for each individual.
#' In the conditional logit model (CLM), the coefficients for a variable are
#' the same for each outcome, because the variables are choice specific.
#' Hence, in the CLM there is a single parameter (a) for each variable z, but
#' there are J values of the variable for each individual.

#' ---------------------------------------------------------------------------
#' ## 6.1 POE5: choice between Pepsi, 7-Up and Coke
#' ---------------------------------------------------------------------------
rm(list=ls())

library(mlogit)

# Cola data
#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/cola.def")

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/cola.rdata"))

# Notice the structure of the data, 3 rows per respondent
head(cola)

# no obs 5466/3=1822 respondents
df <- cola
df$choicename <- rep(c("Pepsi","7-Up","Coke"),1822)

# mean prices and market shares, by brand
df %>% group_by(choicename) %>% summarise(mean.price=mean(price), share=mean(choice))

# Long to wide: one row per respondent
df <- reshape(df, timevar = "choicename", idvar = "id", direction="wide")
head(df)

# The cola dataframe is already stacked, i.e., long.
# We add a variable naming the alternatives, and tell mlogit which
# variable identifies the choice situation (id) and the alternative (brand)
cola$brand <- rep(c("Pepsi","7-Up","Coke"),1822)
DF <- mlogit.data(cola, choice = "choice", shape = "long", chid.var = "id", alt.var = "brand")
head(DF,12)

model1 <- mlogit(choice ~ price, data = DF, reflevel = "Coke") 
summary(model1) # Table 16.4a
head(model.matrix(model1)) # the design matrix

# Predicted probabilities at mean prices
predict(model1)

# Marginal effects at means, formula 16.24 & 16.25
effects(model1, covariate = "price")

# type= first letter is the probability, second is the covariate, a=absolute, r=relative
effects(model1, covariate = "price", type="rr") # elasticity

# Marginal effects at the prices in Table 16.4b
# The rows follow the order of the alternatives in the model: Coke, Pepsi, 7-Up
z <- with(DF, data.frame(price = tapply(price, index(model1)$alt, mean)))
z    # mean prices

z[1,1]=1.1   # Coke
z[2,1]=1     # Pepsi
z[3,1]=1.25  # 7-Up
z

effects(model1, covariate = "price", data = z)
predict(model1, data = z)

#' ### Predicted probabilities and marginal effects with standard errors
b <- as.numeric(coef(model1))
VCVM <- solve(-model1$hessian) # var-cov matrix

names(b) <- c("b11" , "b12", "b2")
b
colnames(VCVM) <- rownames(VCVM) <- names(b)
VCVM

# At given prices
pPepsi <- 1
p7Up <- 1.25
pCoke <- 0.50

# At mean price (overwrites the prices above)
pCoke <- mean(df$price.Coke)
pPepsi <- mean(df$price.Pepsi)
p7Up <- mean(df$`price.7-Up`)

b13 <- 0     # Coke is the reference
b2  <- b[3]
b11 <- b[1]  # pepsi
b12 <- b[2]  # 7Up

# The probability that individual i chooses Pepsi:
PiPepsi <- exp(b11+b2*pPepsi)/
  (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+
     exp(b13+b2*pCoke))
# The probability that individual i chooses Sevenup:
Pi7Up <- exp(b12+b2*p7Up)/
  (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+
     exp(b13+b2*pCoke))
# The probability that individual i chooses Coke:
PiCoke <- 1-PiPepsi-Pi7Up

c(PiCoke,PiPepsi,Pi7Up)

library(CDM)

pred.props <- list("PiPepsi" = ~ I( exp(b11+b2*pPepsi)/
                                      (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))),
                   "Pi7Up" = ~ I(exp(b12+b2*p7Up)/
                                   (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))),
                   "PiCoke" = ~ I(1-(exp(b11+b2*pPepsi)/
                                       (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke)))-
                                    (exp(b12+b2*p7Up)/
                                       (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke)))))

res <- CDM::deltaMethod(pred.props, b, VCVM) # predicted probabilities and se
res

# Marginal effects of the Pepsi price on the three probabilities
dydx_pepsi <- list("dydxPepsi_Pepsi" = ~ I( (exp(b11+b2*pPepsi)/
                                               (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))) *
                                              (1-(exp(b11+b2*pPepsi)/
                                                    (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))))*b2),
                   "dydxPepsi_7Up" = ~ I( (exp(b11+b2*pPepsi)/
                                             (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))) * 
                                            (exp(b12+b2*p7Up)/
                                               (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke)))*b2*-1),
                   "dydxPepsi_Coke" = ~ I( (exp(b11+b2*pPepsi)/
                                              (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke))) *
                                             (1-(exp(b11+b2*pPepsi)/
                                                   (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke)))-
                                                (exp(b12+b2*p7Up)/
                                                   (exp(b11+b2*pPepsi)+exp(b12+b2*p7Up)+exp(b13+b2*pCoke)))) *b2*-1))

resdydx <- CDM::deltaMethod(dydx_pepsi, b, VCVM) # marginal effects and se
resdydx

# Exercise: repeat for
#dydx_7Up
#dydx_cola


#' ---------------------------------------------------------------------------
#' ## 6.2 Cracker data, Franses & Paap
#' ---------------------------------------------------------------------------
rm(list=ls())

library(tidyverse)
library(mlogit)

url <- "https://raw.githubusercontent.com/oysteinm/data/refs/heads/main/cracker.csv"
cracker <- read_csv(url, id=NULL)

# The last purchase of each household is kept out of the estimation
test <- cracker %>% filter(LASTPURCHASE == 0)

# make a choice variable
# PRIVATE = 1
# SUNSHINE = 2
# KEEBLER = 3
# NABISCO = 4

test %>% mutate(choice = 
                  ifelse(PRIVATE==1,"private",
                         ifelse(SUNSHINE==1,"sunshine",
                                ifelse(KEEBLER==1,"keebler","nabisco")))) -> test

test$choice <- as.factor(test$choice)
str(test)

with(test, table(choice,PRIVATE))
with(test, table(choice,SUNSHINE))
with(test, table(choice,KEEBLER))
with(test, table(choice,NABISCO))

names(test) <- c("OBS","HOUSEHOLDID","LASTPURCHASE","PRIVATE","SUNSHINE","KEEBLER","NABISCO",
                 "PRICE.PRIVATE","PRICE.SUNSHINE","PRICE.KEEBLER","PRICE.NABISCO",
                 "DISPL.PRIVATE","DISPL.SUNSHINE","DISPL.KEEBLER","DISPL.NABISCO",
                 "FEAT.PRIVATE","FEAT.SUNSHINE","FEAT.KEEBLER","FEAT.NABISCO",
                 "FEATDISPL.PRIVATE","FEATDISPL.SUNSHINE","FEATDISPL.KEEBLER","FEATDISPL.NABISCO",
                 "choice")

# From wide to long, one variable at a time
test %>% dplyr::select(OBS,HOUSEHOLDID,PRIVATE,SUNSHINE,KEEBLER,NABISCO) -> choice
test %>% dplyr::select(OBS,HOUSEHOLDID,PRICE.PRIVATE,PRICE.SUNSHINE,PRICE.KEEBLER,PRICE.NABISCO) -> price
test %>% dplyr::select(OBS,HOUSEHOLDID,DISPL.PRIVATE,DISPL.SUNSHINE,DISPL.KEEBLER,DISPL.NABISCO) -> displ
test %>% dplyr::select(OBS,HOUSEHOLDID,FEAT.PRIVATE,FEAT.SUNSHINE,FEAT.KEEBLER,FEAT.NABISCO) -> feat
test %>% dplyr::select(OBS,HOUSEHOLDID,FEATDISPL.PRIVATE,FEATDISPL.SUNSHINE,FEATDISPL.KEEBLER,FEATDISPL.NABISCO) -> featdispl

names(choice) <- c("obs","householdid","private","sunshine","keebler","nabisco")
choice %>% gather(brand, choice, private, sunshine, keebler, nabisco) %>% arrange(obs) -> choice

names(price) <- c("obs","householdid","private","sunshine","keebler","nabisco")
price %>% gather(brand, price, private, sunshine, keebler, nabisco) %>% arrange(obs) -> price

names(displ) <- c("obs","householdid","private","sunshine","keebler","nabisco")
displ %>% gather(brand, displ, private, sunshine, keebler, nabisco) %>% arrange(obs) -> displ

names(feat) <- c("obs","householdid","private","sunshine","keebler","nabisco")
feat %>% gather(brand, feat, private, sunshine, keebler, nabisco) %>% arrange(obs) -> feat

names(featdispl) <- c("obs","householdid","private","sunshine","keebler","nabisco")
featdispl %>% gather(brand, featdispl, private, sunshine, keebler, nabisco) %>% arrange(obs) -> featdispl

temp1 <- merge(choice,price, by=c("obs","householdid","brand"))
temp2 <- merge(temp1,displ, by=c("obs","householdid","brand"))
temp3 <- merge(temp2,feat, by=c("obs","householdid","brand"))
df.cracker <- merge(temp3,featdispl, by=c("obs","householdid","brand")) 
df.cracker %>% arrange(obs) -> df.cracker

rm(temp1,temp2,temp3,test,displ,feat,featdispl,price,choice)

df.cracker %>% group_by(brand) %>% summarise(mean.price=mean(price), share=mean(choice))

df.cracker <- mlogit.data(df.cracker, choice = "choice", shape = "long",
                          chid.var = "obs", alt.var = "brand")

names(df.cracker)
head(df.cracker)

# conditional logit model
model1 <- mlogit(choice ~ displ+feat+featdispl+price, data = df.cracker, reflevel = "nabisco") 
summary(model1) # p. 101 Franses & Paap
head(model.matrix(model1)) # the design matrix

# Predicted values at means
predict(model1)

# Marginal effects at means
effects(model1, covariate = "price")

# first letter in type= probability, second is covariate, a=absolute, r=relative
effects(model1, covariate = "price", type="rr") # elasticity


#' ---------------------------------------------------------------------------
#' ## 6.3 Nested logit model, Cracker data
#' ---------------------------------------------------------------------------
#' The private label is in its own nest, the three national brands in another.
model2 <- mlogit(choice ~ displ+feat+featdispl+price, df.cracker, reflevel = "nabisco",
                 nests = list(private="private", brand=c("keebler","nabisco","sunshine")),
                 unscaled=TRUE) 
summary(model2) 

# Compare the conditional and nested logit, likelihood ratio test
lrtest(model1, model2)

# NB: predict() and effects() did not work for the nested model:
# predict(model2)
# effects(model2, covariate = "price")
# effects(model2, covariate = "price", type="rr") # elasticity


#' ---------------------------------------------------------------------------
#' ## 6.4 Mixing individual- and alternative-specific variables (Fishing data)
#' ---------------------------------------------------------------------------
#browseURL("https://cran.r-project.org/web/packages/mlogit/index.html")
#' The mlogit formula has three parts:  y ~ a | b | c
#'   a: alternative-specific variables with a generic coefficient (as in the CLM)
#'   b: individual-specific variables (as in the MNLM)
#'   c: alternative-specific variables with alternative-specific coefficients
rm(list=ls())
library(mlogit)

data("Fishing", package = "mlogit")
head(Fishing)

Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
m <- mlogit(mode ~ price | income | catch, data = Fish)
summary(m)

# compute a data.frame containing the mean value of the covariates in the sample
z <- with(Fish, data.frame(price = tapply(price, index(m)$alt, mean), 
                           catch = tapply(catch, index(m)$alt, mean), 
                           income = mean(income)))
z

# compute the marginal effects (the second one is an elasticity, see the
# effects.mlogit help page: https://cran.r-project.org/web/packages/mlogit/mlogit.pdf)
effects(m, covariate = "income", data = z)
# elasticity
effects(m, covariate = "price", type = "rr", data = z)


#' ===========================================================================
#' # Part 7: Ordered probit and logit models
#' ===========================================================================
rm(list=ls())

library(MASS)
library(tidyverse)

#' Data definition file:
#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/nels_small.def")
#' Read the data
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/nels_small.rdata"))

# psechoice:
# = 1 if first postsecondary education was no college
# = 2 if first postsecondary education was a 2-year college
# = 3 if first postsecondary education was a 4-year college
# 
# grades		= average grade in math, english and social studies on 13 point scale with 1 = highest

#' Here we treat psechoice as ordered: more education is "more".

# Descriptive statistics
table(nels_small$psechoice)
histogram(~psechoice, data = nels_small, type="percent", label=TRUE)
favstats(~psechoice, data = nels_small)

favstats(~grades, data = nels_small)
densityplot(~grades, data = nels_small)

# The dependent variable must be a factor
nels_small$psechoiceF <- factor(nels_small$psechoice, labels=c("no college","2-year college","4-year college"))

#' Ordered logit model
# polr stands for proportional odds logistic regression (ordered logit)
ologit <- polr(psechoiceF ~ grades, data=nels_small, Hess=TRUE, method = "logistic")
summary(ologit)
# In one sweep
summary(polr(as.factor(psechoice) ~ grades, data=nels_small, Hess=TRUE, method = "logistic"))

#' Ordered probit
oprobit <- polr(psechoiceF ~ grades, data=nels_small, Hess=TRUE, method="probit")
summary(oprobit)

## store table
ctable <- coef(summary(ologit))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
ctable <- cbind(ctable, "p-value" = p)
ctable

confint(ologit)
confint.default(ologit) # CIs assuming normality

# Ordered logit model odds ratio
exp(coef(ologit))

library(stargazer)
stargazer(ologit, oprobit, type="text")

#' Ordered logit predicted probabilities
# At the mean of grades
xmeans <- nels_small %>% dplyr::select(grades) %>% summarise(grades=mean(grades))
fitted <- predict(ologit, newdata=xmeans, type="probs")
fitted

# For all values of grades
newdat <- data.frame(grades = 1:13)
newdat <- cbind(newdat, predict(ologit, newdata=newdat, type = "probs"))
head(newdat)

# Wide to long
lnewdat <- newdat %>% gather(Level, Probability, -grades)
## view first few rows
head(lnewdat)

ggplot(lnewdat, aes(x = grades, y = Probability, colour = Level)) +  geom_line()

# First differences: change in the probabilities when grades go from 3 to 4
pr34 <- predict(ologit, newdata=data.frame(grades=c(3,4)), type="probs")
pr34
pr34[2,]-pr34[1,]

#' Marginal effects
#browseURL("https://cran.r-project.org/src/contrib/Archive/oglmx/")
# Install from .tar.gz locally
library(oglmx)

results.oprob <- oprobit.reg(psechoice ~ grades, data=nels_small)
summary(results.oprob)

margins.oglmx(results.oprob, atmeans = TRUE, AME = FALSE, location = NULL, outcomes = "All",
              ascontinuous = FALSE, Vars = NULL)

#' ### More marginal effects, the housing data (MASS)
library(erer)
data(housing) # press F1 for info on data
str(housing); head(housing)

# Fit an ordered choice model with polr from the MASS library.
# The original specification used in MASS is Sat ~ Infl + Type + Cont,
# Freq is added here as a continuous variable so we can plot
# the probabilities against it.
fm <- Sat ~ Infl + Type + Cont + Freq
ra <- polr(fm, data = housing, weights = Freq, Hess = TRUE, method = "probit")
rb <- polr(fm, data = housing, weights = Freq, Hess = TRUE, method = "logistic")
summary(ra); summary(rb)

# Compute the marginal effect
mea <- ocME(w = ra); mea
meb <- ocME(w = rb); meb
mea$out
meb$out

# Predicted probabilities as a function of Freq
fa <- ocProb(w = ra, nam.c = "Freq", n = 300)
fb <- ocProb(w = rb, nam.c = "Freq", n = 300)
plot(fa)
plot(fb) 

# slightly different plot
library(effects)
plot(effect("Freq", rb))


#' ===========================================================================
#' # Part 8: Censored data, the Tobit model
#' ===========================================================================
#' Example 16.16
rm(list=ls())

# https://stats.idre.ucla.edu/r/dae/tobit-models/
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/mroz.rdata"))

head(mroz)

# hours worked is censored at zero
library(VGAM)

summary(fit <- vglm(hours ~ educ + exper + age + kidsl6, tobit(Lower = 0), data = mroz))

library(censReg)

summary(fit1 <- censReg(hours ~ educ + exper + age + kidsl6, left = 0, data = mroz))

# Marginal effects at the means, with kidsl6 = 1
margEff(fit1, xValues = c(1, mean(mroz$educ), mean(mroz$exper), mean(mroz$age), 1))

summary(margEff(fit1, xValues = c(1, mean(mroz$educ), mean(mroz$exper), mean(mroz$age), 1)))

# manually marginal effect of educ: b*Phi(x'b/sigma)
mean.mroz <- data.frame(
  educ = mean(mroz$educ),
  exper = mean(mroz$exper),
  age = mean(mroz$age),
  kidsl6 = 1
)

p <- predict(fit, newdata = mean.mroz)

fit1$estimate[2]*pnorm(p[1]/exp(fit1$estimate[6])) # estimate[6] is log(sigma)

#' ### Right censoring, simulated data
library(rockchalk)

set.seed(421)
df1 <- genCorrelatedData2(N = 159, means = c(6.02, 5.07, 5.12), sds = c(0.4, 0.4, 0.5),
                          rho = c(0,0,0), stde = 0.5, beta = c(0, 0.992, 0.616, -0.292))

summary(lm(y~x1+x2+x3, data=df1))

library(tidyverse)

# censor y from above at 10
df1 <- df1 %>% 
  mutate(y_right=ifelse(y>9,10,y))

summary(fit2 <- censReg(y_right ~ x1 + x2 + x3, right = 10, data = df1))

margEff(fit2)

summary(margEff(fit2))


#' ===========================================================================
#' # Part 9: Extension, a Bayesian logit model (coke data, Example 16.6)
#' ===========================================================================
rm(list=ls())

library(tidyverse)

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/coke.rdata"))

# The frequentist model, for comparison
logit <- glm(coke~pratio+disp_coke+disp_pepsi, family = binomial(link = "logit"), data = coke)
summary(logit)

library(brms)

#' Priors for the beta's: a Student-t with 3 df, location 0 and scale 2.5,
#' a weakly informative default prior for logistic regression, see
#' Gelman, Jakulin, Pittau & Su (2008), Annals of Applied Statistics
#browseURL("https://projecteuclid.org/download/pdfview_1/euclid.aoas/1231424214")

x_values <- seq(-15,15, length.out = 1000)

data.frame(x_values) %>%
  ggplot(aes(x_values))+
  stat_function(fun = function(x) dt(x/2.5, df=3)/2.5)

# brms
eqtn <- bf(coke~pratio+disp_pepsi+disp_coke)

prior1 <- get_prior(eqtn, family = "bernoulli", data=coke)
prior1
prior1$prior[c(2,3,4)] <- "student_t(3, 0, 2.5)"
prior1

# This takes a few minutes to run, depending on your computer!
fit_brms <- brm(eqtn, data = coke,
                cores=4, chains = 4, iter=1000, warmup=200,
                family = "bernoulli",
                prior = prior1)

summary(fit_brms)

plot(fit_brms)
plot(conditional_effects(fit_brms), points = TRUE)

# Posterior summaries: mean, sd and 95% credible intervals
posterior_summary(fit_brms)

# Posterior predictive check
pp <- pp_check(fit_brms)
pp # similar density plot between observed and predicted

## perform two-sided hypothesis testing
hyp1 <- hypothesis(fit_brms, "disp_pepsi=disp_coke")
hyp1
plot(hyp1)

#' Exercise: compare the posterior means with the ML estimates from glm().
#' Why are they so similar?
