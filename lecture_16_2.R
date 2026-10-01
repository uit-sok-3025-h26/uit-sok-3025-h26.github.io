# Lecture 16.2 — original nels_small and cola data
# Packages: mlogit, nnet, MASS, marginaleffects (dfidx is an mlogit dependency)
# Run interactively, or Rscript lecture_16_2.R from this folder.
dir.create("figures", showWarnings = FALSE)

# ---- setup ----
library(mlogit)
library(nnet)
library(marginaleffects)

# Definitions:
# https://www.principlesofeconometrics.com/poe5/data/def/nels_small.def
# https://www.principlesofeconometrics.com/poe5/data/def/cola.def
load(url("https://www.principlesofeconometrics.com/poe5/data/rdata/nels_small.rdata"))
load(url("https://www.principlesofeconometrics.com/poe5/data/rdata/cola.rdata"))
head(nels_small)
head(cola)

names_education <- c("No college", "Two-year", "Four-year")
nels_small$education <- factor(nels_small$psechoice, levels = 1:3,
                               labels = names_education)
nels_small$education_ordered <- ordered(nels_small$education,
                                        levels = names_education)
# Larger grades means WORSE performance. No college is the reference.
# multinom estimates the same MNL as mlogit(education ~ 0 | grades).
# This fitted class supports marginaleffects::slopes() directly.
mnl <- nnet::multinom(education ~ grades, data = nels_small,
                      trace = FALSE, Hess = TRUE, maxit = 1000,
                      reltol = 1e-10)
summary(mnl)
prop.table(table(nels_small$education))
colMeans(fitted(mnl))

# Cola records are in Pepsi, 7-Up, Coke order WITHIN each original id.
# Assign brand labels BEFORE any sorting; id identifies a complete choice set.
names_drinks <- c("Pepsi", "7-Up", "Coke")
cola$brand <- factor(rep(names_drinks, times = nrow(cola) / 3),
                     levels = names_drinks)
stopifnot(all(table(cola$id) == 3),
          all(tapply(cola$choice, cola$id, sum) == 1))
cola_long <- mlogit.data(cola, choice = "choice", shape = "long",
                         chid.var = "id", alt.var = "brand")
cl <- mlogit(choice ~ price | 1, data = cola_long, reflevel = "Coke")
summary(cl)
# The price-only specification matches the main textbook example.
# feature and display are retained in cola for an optional extension.

# Stable softmax and log-sum; all coefficients below come from fitted models.
softmax <- function(v) {
  w <- exp(v - max(v))
  w / sum(w)
}
logsum <- function(v) max(v) + log(sum(exp(v - max(v))))
B <- coef(mnl)
alpha_edu <- c(0, B["Two-year", "(Intercept)"], B["Four-year", "(Intercept)"])
beta_edu <- c(0, B["Two-year", "grades"], B["Four-year", "grades"])
education_prob <- function(grades) {
  setNames(softmax(alpha_edu + beta_edu * grades), names_education)
}
education_me <- function(grades) {
  p <- education_prob(grades)
  p * (beta_edu - sum(p * beta_edu))
}
b_cl <- coef(cl)
alpha_drink <- c(unname(b_cl["(Intercept):Pepsi"]),
                 unname(b_cl["(Intercept):7-Up"]), 0)
beta_price <- unname(b_cl["price"])
stopifnot(!anyNA(c(alpha_drink, beta_price)))
drink_index <- function(prices) alpha_drink + beta_price * prices
drink_prob <- function(prices) {
  setNames(softmax(drink_index(prices)), names_drinks)
}

# ---- by-hand ----
v <- alpha_edu + beta_edu * 6.64
rbind(original = softmax(v),               # R1: location invariance
      shifted  = softmax(v + 7),
      scaled   = softmax(2 * v))           # scale is NOT innocuous: beta/sigma

h <- 1e-6                                  # R2: numerical gradient and Hessian
grad_L <- sapply(1:3, function(j) (logsum(v + h * (1:3 == j)) - logsum(v)) / h)
p <- softmax(v)
Sigma <- diag(p) - p %*% t(p)
list(gradient_vs_p = rbind(grad_L, p), Sigma = Sigma, row_sums = rowSums(Sigma))

c(two_vs_none = beta_edu[2], four_vs_none = beta_edu[3],   # R3: every pair
  two_vs_four = beta_edu[2] - beta_edu[3])
exp(beta_edu[2] - beta_edu[3])   # odds of two-year vs four-year, per GRADES unit

# ---- education ----
education_prob(2.635)  # better grades
education_prob(6.64)   # median grades
rbind(`ME at 2.635` = education_me(2.635), `ME at 6.64` = education_me(6.64))
g_star <- optimize(function(x) education_prob(x)[2],
                   interval = range(nels_small$grades), maximum = TRUE)$maximum
# An interior maximum has a zero slope; a boundary maximum need not.
c(g_star = g_star, education_prob(g_star))
g <- seq(min(nels_small$grades), max(nels_small$grades), length.out = 200)
edu_matrix <- t(vapply(g, education_prob, numeric(3)))
matplot(g, edu_matrix, type = "l", lty = 1, lwd = 2,
        col = c("#555555", "#2487a8", "#d45d42"), ylim = c(0, 1),
        xlab = "GRADES (larger = worse performance)",
        ylab = "Predicted choice probability")
abline(v = g_star, lty = 2, col = "#2487a8")
legend("right", names_education, lty = 1, lwd = 2,
       col = c("#555555", "#2487a8", "#d45d42"), bty = "n")
png("figures/education.png", width = 1000, height = 620, res = 140)
matplot(g, edu_matrix, type = "l", lty = 1, lwd = 2,
        col = c("#555555", "#2487a8", "#d45d42"), ylim = c(0, 1),
        xlab = "GRADES (larger = worse performance)",
        ylab = "Predicted choice probability")
abline(v = g_star, lty = 2, col = "#2487a8")
legend("right", names_education, lty = 1, lwd = 2,
       col = c("#555555", "#2487a8", "#d45d42"), bty = "n")
dev.off()

# ---- education-elasticity ----
grade_profiles <- data.frame(grades = c(2.635, 6.64))
slopes(mnl, variables = "grades", newdata = grade_profiles, slope = "eyex")
avg_slopes(mnl, variables = "grades", slope = "eyex")
rbind(at_2.635 = education_me(2.635) * 2.635 / education_prob(2.635),
      at_6.64 = education_me(6.64) * 6.64 / education_prob(6.64))

# ---- education-differences ----
change <- education_prob(6.64) - education_prob(2.635)
change
sum(change)              # zero up to floating point error: finite changes
sum(education_me(6.64))  # zero: derivatives (rows of Sigma sum to zero)

# ---- drink-own-price ----
base_prices <- c(Pepsi = 1.00, `7-Up` = 1.25, Coke = 1.10)
pepsi_up <- base_prices
pepsi_up["Pepsi"] <- 1.10
drink_table <- rbind(Baseline = drink_prob(base_prices),
                     `Pepsi +10 cents` = drink_prob(pepsi_up))
round(drink_table, 4)
round(drink_table[2, ] - drink_table[1, ], 4)
barplot(t(drink_table), beside = TRUE, ylim = c(0, 0.6),
        col = c("#2876a8", "#68a56c", "#c85e58"),
        ylab = "Choice probability", legend.text = names_drinks,
        args.legend = list(x = "topright", bty = "n"))
png("figures/drink-own-price.png", width = 1000, height = 620, res = 140)
barplot(t(drink_table), beside = TRUE, ylim = c(0, 0.6),
        col = c("#2876a8", "#68a56c", "#c85e58"),
        ylab = "Choice probability", legend.text = names_drinks,
        args.legend = list(x = "topright", bty = "n"))
dev.off()

p0 <- drink_prob(base_prices)
jacobian <- beta_price * (diag(p0) - p0 %*% t(p0))           # R5
dimnames(jacobian) <- list(prob = names_drinks, price = names_drinks)
round(jacobian, 4)       # first row reproduces Table 16.4b
elasticity <- jacobian * outer(1 / p0, base_prices)          # dp_j/dP_k * P_k/p_j
round(elasticity, 3)     # off-diagonal columns are constant: IIA

# ---- cola-elasticity-package ----
# A symmetric log-price perturbation estimates d log(p_j) / d log(P_Pepsi).
# id1/id2 are the choice-set and alternative indices expected by the
# marginaleffects mlogit adapter (see the linked package vignette).
h_log <- 1e-4
nd <- as.data.frame(cola_long)
idx <- dfidx::idx(cola_long)
nd$id1 <- idx[[1]]
nd$id2 <- idx[[2]]
nd$idx <- NULL
# Give every shopper the same stated price profile for this comparison.
nd$price <- unname(base_prices[as.character(nd$id2)])
# The marginaleffects mlogit adapter reads rows in the model's internal
# alternative order (reference level first: Coke, Pepsi, 7-Up), not by id2.
# Without this sort, prices are silently attached to the wrong brands.
alt_order <- names(cl$freq)
nd <- nd[order(nd$id1, match(as.character(nd$id2), alt_order)), ]
# Guard: baseline predictions must reproduce drink_prob(base_prices).
base_check <- aggregate(estimate ~ group, FUN = mean,
                        data = predictions(cl, newdata = nd))
stopifnot(isTRUE(all.equal(
  base_check$estimate,
  unname(drink_prob(base_prices)[as.character(base_check$group)]),
  tolerance = 1e-6)))
lo <- transform(nd, term = "lo",
                price = ifelse(id2 == "Pepsi", price * exp(-h_log), price))
hi <- transform(nd, term = "hi",
                price = ifelse(id2 == "Pepsi", price * exp(h_log), price))
scenarios <- rbind(lo, hi)
scenarios$term <- factor(scenarios$term, levels = c("lo", "hi"))

# All shoppers have the same profile here, so group means equal each
# profile probability. Explicit matching avoids depending on row order.
pepsi_elasticity <- function(x) {
  a <- aggregate(estimate ~ term + group, data = x, FUN = mean)
  low <- a[a$term == "lo", ]
  high <- a[a$term == "hi", ]
  high <- high[match(low$group, high$group), ]
  data.frame(term = as.character(low$group),
             estimate = (log(high$estimate) - log(low$estimate)) / (2 * h_log))
}
predictions(cl, newdata = scenarios, hypothesis = pepsi_elasticity)
# Compare the three estimates with the Pepsi-price column of the exact matrix.
elasticity[, "Pepsi"]

# ---- drink-cross-price ----
coke_up <- base_prices
coke_up["Coke"] <- 1.25
round(drink_prob(coke_up) - drink_prob(base_prices), 4)

# R6: after/before ratios for a Pepsi price increase
round(drink_prob(pepsi_up) / drink_prob(base_prices), 4)
exp(logsum(drink_index(base_prices)) - logsum(drink_index(pepsi_up)))  # S_old / S_new

# ---- iia ----
sevenup_up <- base_prices
sevenup_up["7-Up"] <- 1.40
before <- drink_prob(base_prices)
after <- drink_prob(sevenup_up)
rbind(before = before, after = after)
c(before_ratio = before["Pepsi"] / before["Coke"],
  after_ratio = after["Pepsi"] / after["Coke"])

# Red bus / blue bus with a nested logit (Appendix B): car alone; buses in one nest.
nested_bus <- function(lambda) {
  bus_nest <- 2^lambda                    # (e^0 + e^0)^lambda
  c(car = 1, red_bus = bus_nest / 2, blue_bus = bus_nest / 2) / (1 + bus_nest)
}
round(sapply(c(`MNL (lambda = 1)` = 1, `lambda = 0.5` = 0.5,
               `lambda = 0.01` = 0.01), nested_bus), 3)

# ---- welfare ----
alpha_money <- -beta_price                 # marginal utility of one dollar
c(delta_CS  = (logsum(drink_index(pepsi_up)) -
               logsum(drink_index(base_prices))) / alpha_money,
  first_order = -drink_prob(base_prices)[["Pepsi"]] * 0.10)

# ---- ordered ----
op <- MASS::polr(education_ordered ~ grades, data = nels_small,
                 method = "probit", Hess = TRUE)
summary(op)
beta_op <- unname(coef(op)["grades"])
mu1 <- unname(op$zeta[1]); mu2 <- unname(op$zeta[2])
ordered_prob <- function(grades) {
  low <- pnorm(mu1 - beta_op * grades)
  below_high <- pnorm(mu2 - beta_op * grades)
  setNames(c(low, below_high - low, 1 - below_high), names_education)
}
ordered_me <- function(grades) {                                   # R8
  a1 <- mu1 - beta_op * grades; a2 <- mu2 - beta_op * grades
  setNames(beta_op * c(-dnorm(a1), dnorm(a1) - dnorm(a2), dnorm(a2)),
           names_education)
}
x_star <- (mu1 + mu2) / (2 * beta_op)
ordered_matrix <- t(vapply(g, ordered_prob, numeric(3)))
matplot(g, ordered_matrix, type = "l", lty = 1, lwd = 2,
        col = c("#555555", "#2487a8", "#d45d42"), ylim = c(0, 1),
        xlab = "GRADES (larger = worse performance)",
        ylab = "Ordered-probit probability")
abline(v = x_star, lty = 2, col = "#2487a8")
legend("right", names_education, lty = 1, lwd = 2,
       col = c("#555555", "#2487a8", "#d45d42"), bty = "n")
png("figures/ordered.png", width = 1000, height = 620, res = 140)
matplot(g, ordered_matrix, type = "l", lty = 1, lwd = 2,
        col = c("#555555", "#2487a8", "#d45d42"), ylim = c(0, 1),
        xlab = "GRADES (larger = worse performance)",
        ylab = "Ordered-probit probability")
abline(v = x_star, lty = 2, col = "#2487a8")
legend("right", names_education, lty = 1, lwd = 2,
       col = c("#555555", "#2487a8", "#d45d42"), bty = "n")
dev.off()
ordered_prob(2.635)
ordered_prob(6.64)
rbind(`OP ME at 2.635` = ordered_me(2.635), `OP ME at 6.64` = ordered_me(6.64))
c(x_star_ordered = x_star, g_star_mnl = g_star)

# R8 (ii): slopes of qnorm(P(Y <= k)) in GRADES. Ordered probit: -beta for both k.
cum_slope <- function(prob_fun, x, h = 1e-5) {
  cz <- function(x) { p <- prob_fun(x); qnorm(c(p[1], p[1] + p[2])) }
  setNames((cz(x + h) - cz(x)) / h, c("P(Y<=1)", "P(Y<=2)"))
}
rbind(ordered_2.635 = cum_slope(ordered_prob, 2.635),
      ordered_6.64  = cum_slope(ordered_prob, 6.64),
      mnl_2.635     = cum_slope(education_prob, 2.635),
      mnl_6.64      = cum_slope(education_prob, 6.64))

# ---- ordered-elasticity ----
slopes(op, variables = "grades", newdata = grade_profiles, slope = "eyex")
avg_slopes(op, variables = "grades", slope = "eyex")
rbind(at_2.635 = ordered_me(2.635) * 2.635 / ordered_prob(2.635),
      at_6.64 = ordered_me(6.64) * 6.64 / ordered_prob(6.64))

# ---- estimation-check ----
# Rows are choice sets, columns are explicitly ordered by brand.
P <- fitted(cl, outcome = FALSE)[, names_drinks, drop = FALSE]
Y <- xtabs(choice ~ id + brand, data = cola)[, names_drinks]
price_matrix <- xtabs(price ~ id + brand, data = cola)[, names_drinks]
# Align with fitted choice-set row names rather than relying on sort order.
Y <- Y[rownames(P), , drop = FALSE]
price_matrix <- price_matrix[rownames(P), , drop = FALSE]
rbind(observed_shares = colMeans(Y), fitted_shares = colMeans(P))
c(observed_price_paid = mean(rowSums(Y * price_matrix)),
  fitted_price_paid = mean(rowSums(P * price_matrix)))
vcov(cl) # Full estimation covariance matrix, available from actual fitting.
