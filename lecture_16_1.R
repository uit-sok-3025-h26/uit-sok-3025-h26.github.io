# Lecture 16.1 — Binary choices using the original Coke data
# Requires: tidyverse, marginaleffects, this.path
# Companion script to lecture_16_1.qmd: same code, plus figures saved to figures/.

# Install if you don't have it: install.packages("this.path")
library(this.path)

# Sets directory to the folder containing the script
setwd(this.dir())
getwd()

# Run from this folder. Figures are saved in figures/.
dir.create("figures", showWarnings = FALSE)

# ---- setup ----
library(tidyverse)
library(marginaleffects)

# Original textbook observations: 1 = Coke chosen, 0 = Pepsi chosen.
# Definition: https://www.principlesofeconometrics.com/poe5/data/def/coke.def
# Download once and reuse the local copy, so rendering does not depend on the server.
data_file <- "coke.rdata"
if (!file.exists(data_file)) {
  download.file("https://www.principlesofeconometrics.com/poe5/data/rdata/coke.rdata",
                data_file, mode = "wb")
}
load(data_file)
coke <- coke %>% mutate(pratio = pr_coke / pr_pepsi)
head(coke)
summary(coke)

# Plot only where the data are: 1st to 99th percentile of the price ratio.
price_range <- quantile(coke$pratio, c(0.01, 0.99))
coke_in_range <- filter(coke, between(pratio, price_range[[1]], price_range[[2]]))
price_range

# pratio: Coke price / Pepsi price; displays are binary indicators.
# Keep this regressor order for the by-hand coefficient calculations.
probit <- glm(coke ~ pratio + disp_coke + disp_pepsi,
              family = binomial(link = "probit"), data = coke)
logit <- glm(coke ~ pratio + disp_coke + disp_pepsi,
             family = binomial(link = "logit"), data = coke)
summary(logit)
summary(probit)

# ---- pdf-cdf-figure ----
z <- seq(-3.5, 3.5, length.out = 300)
curves <- rbind(
  data.frame(z, value = dnorm(z), function_name = "Normal PDF"),
  data.frame(z, value = pnorm(z), function_name = "Normal CDF"),
  data.frame(z, value = dlogis(z), function_name = "Logistic PDF"),
  data.frame(z, value = plogis(z), function_name = "Logistic CDF")
)
figure_plot <- ggplot(curves, aes(z, value, colour = function_name)) +
  geom_line(linewidth = 1) + facet_wrap(~ function_name, scales = "free_y") +
  labs(x = "Choice index / threshold", y = "Density or cumulative probability") +
  theme_minimal()
print(figure_plot)
ggsave("figures/pdf-cdf-figure.png", figure_plot, width = 8, height = 4.5, dpi = 150)

# ---- by-hand ----
b <- coef(logit)
coef(logit) / coef(probit)                         # R1: scale normalization

c(logit  = -b[["disp_coke"]] / b[["pratio"]],       # R1: identified MRS
  probit = -coef(probit)[["disp_coke"]] / coef(probit)[["pratio"]])
hypotheses(logit, hypothesis = "-disp_coke / pratio = 0")   # with a delta-method SE

exp(b)                                              # R2: odds ratios

# The reference shopper used throughout: ratio 1.1, no displays.
profile0 <- data.frame(pratio = 1.1, disp_coke = 0, disp_pepsi = 0)
x0 <- c(1, profile0$pratio, profile0$disp_coke, profile0$disp_pepsi)
eta0 <- sum(x0 * b); p0 <- plogis(eta0)
c(eta = eta0, p = p0,
  slope_by_hand = p0 * (1 - p0) * b[["pratio"]],        # R3
  slope_bound   = -0.25 * abs(b[["pratio"]]),           # R4
  elasticity    = (1 - p0) * b[["pratio"]] * x0[2])     # R6

# ---- choice-curves ----
grid <- data.frame(pratio = seq(price_range[[1]], price_range[[2]], length.out = 120),
                   disp_coke = 0, disp_pepsi = 0)
plotdat <- bind_rows(
  Logit  = as_tibble(predictions(logit,  newdata = grid)),
  Probit = as_tibble(predictions(probit, newdata = grid)),
  .id = "model")
figure_plot <- ggplot(plotdat, aes(pratio, estimate, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = coke_in_range, aes(x = pratio), inherit.aes = FALSE,
           sides = "b", alpha = 0.1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Coke price / Pepsi price", y = "P(buy Coke)",
       colour = NULL, fill = NULL) + theme_minimal()
print(figure_plot)
ggsave("figures/choice-curves.png", figure_plot, width = 8, height = 4.5, dpi = 150)

# The reference shopper (ratio 1.1, no displays) under both links.
bind_rows(Logit  = as_tibble(predictions(logit,  newdata = profile0)),
          Probit = as_tibble(predictions(probit, newdata = profile0)),
          .id = "model") %>%
  select(model, estimate, conf.low, conf.high)   # interval is computed on the link scale
                                                # and back-transformed, so no SE column

# ---- price-slopes ----
profiles <- data.frame(pratio = c(0.8, 1.1, 1.5),
                       disp_coke = 0, disp_pepsi = 0)
slopes(logit, variables = "pratio", newdata = profiles)
avg_slopes(logit, variables = "pratio")                  # AME
slopes(logit, variables = "pratio", newdata = "mean")    # MEM

# R7 by hand: AME averages f(eta_i); MEM evaluates f at the mean index.
p_i <- fitted(logit)
x_bar <- c(1, colMeans(coke[, c("pratio", "disp_coke", "disp_pepsi")]))
p_bar <- plogis(sum(x_bar * b))
c(AME = mean(p_i * (1 - p_i)) * b[["pratio"]],
  MEM = p_bar * (1 - p_bar) * b[["pratio"]])

# Is the whole index inside the region where the logistic density is concave?
eta_i <- predict(logit, type = "link")
c(min_eta = min(eta_i), max_eta = max(eta_i),
  concave_bound = log(2 + sqrt(3)),
  share_in_concave_region = mean(abs(eta_i) < log(2 + sqrt(3))))

# ---- finite-price-change ----
price_pair <- data.frame(pratio = c(1.0, 1.1), disp_coke = 0, disp_pepsi = 0)
p_pair <- predict(logit, price_pair, type = "response")
data.frame(price_pair, probability = p_pair)

# R5: exact change versus first- and second-order approximations at x = 1.0
p1 <- p_pair[[1]]; dx <- 0.1; bP <- b[["pratio"]]
c(exact        = unname(p_pair[2] - p_pair[1]),
  tangent      = p1 * (1 - p1) * bP * dx,
  second_order = p1 * (1 - p1) * bP * dx +
                 0.5 * p1 * (1 - p1) * (1 - 2 * p1) * (bP * dx)^2)

# ---- display-comparison ----
display_pair <- data.frame(pratio = 1.1, disp_coke = c(0, 1), disp_pepsi = 0)
predict(logit, display_pair, type = "response")
diff(predict(logit, display_pair, type = "response"))
avg_comparisons(logit, variables = "disp_coke")

# ---- elasticity ----
scales <- c("dydx", "dyex", "eydx", "eyex")

# All four scales at stated profiles, for both links, in one table.
elas <- imap_dfr(list(Logit = logit, Probit = probit), \(fit, model_name)
  map_dfr(scales, \(scale_name)
    slopes(fit, variables = "pratio", type = "response",
           slope = scale_name, newdata = profiles) %>%
      as_tibble() %>%
      transmute(model = model_name, scale = scale_name, pratio, estimate)))
elas %>%
  pivot_wider(names_from = scale, values_from = estimate) %>%
  mutate(across(all_of(scales), \(v) round(v, 3)))

# Average observation-specific effects over the original sample (logit).
map_dfr(scales, \(scale_name)
  avg_slopes(logit, variables = "pratio", type = "response",
             slope = scale_name) %>%
    as_tibble() %>%
    transmute(scale = scale_name, estimate, conf.low, conf.high))

# Verify all four logit measures by hand at the same profiles.
p_profile <- predict(logit, profiles, type = "response")
m_profile <- p_profile * (1 - p_profile) * coef(logit)[["pratio"]]
data.frame(profiles, probability = p_profile,
           dydx = m_profile, dyex = m_profile * profiles$pratio,
           eydx = m_profile / p_profile,
           eyex = m_profile * profiles$pratio / p_profile)

# Exact 1% price-ratio increase versus the local elasticity approximation.
p_plus <- predict(logit, transform(profiles, pratio = 1.01 * pratio),
                  type = "response")
data.frame(profiles,
           exact_percent = 100 * (p_plus / p_profile - 1),
           local_percent = m_profile * profiles$pratio / p_profile,
           exact_percentage_points = 100 * (p_plus - p_profile))

# ---- elasticity-curves ----
# Price elasticity across all four display combinations.
elasticity_grid <- expand.grid(
  pratio = seq(price_range[[1]], price_range[[2]], length.out = 150),
  disp_coke = 0:1, disp_pepsi = 0:1)
elasticity_grid$elasticity <- slopes(
  logit, variables = "pratio", slope = "eyex", type = "response",
  newdata = elasticity_grid)$estimate
elasticity_grid$display <- with(elasticity_grid,
  paste0("Coke = ", disp_coke, ", Pepsi = ", disp_pepsi))
figure_plot <- ggplot(elasticity_grid, aes(pratio, elasticity, colour = display)) +
  geom_line(linewidth = 1) +
  labs(x = "Coke price / Pepsi price", y = "Choice-probability elasticity",
       colour = "Displays") + theme_minimal()
print(figure_plot)
ggsave("figures/elasticity-curves.png", figure_plot, width = 8, height = 4.5, dpi = 150)

# ---- induced-interaction ----
slope_grid <- expand.grid(pratio = c(0.9, 1.3),
                          disp_coke = 0:1, disp_pepsi = 0)
# Additive model: no product term, yet slopes differ by display status.
cbind(slope_grid,
      slope = slopes(logit, variables = "pratio", newdata = slope_grid)$estimate,
      eta0  = b[["(Intercept)"]] + b[["pratio"]] * slope_grid$pratio)
c(threshold = -b[["disp_coke"]] / 2)   # R8: sign flips where eta0 = -beta_C / 2

# The price ratio at which eta0 = -beta_C / 2 (no displays for Pepsi).
c(pratio_flip = (-b[["disp_coke"]] / 2 - b[["(Intercept)"]]) / b[["pratio"]])

# ---- interaction ----
# Product term with the price ratio centred at parity (pratio = 1).
logit_interaction <- glm(coke ~ I(pratio - 1) * disp_coke + disp_pepsi,
                         family = binomial(link = "logit"), data = coke)
logit_interaction_raw <- glm(coke ~ pratio * disp_coke + disp_pepsi,
                             family = binomial(link = "logit"), data = coke)

# Centring changes only the intercept and the display main effect.
coef_names <- c("(Intercept)", "pratio", "disp_coke", "disp_pepsi", "pratio:disp_coke")
round(rbind(uncentred = setNames(coef(logit_interaction_raw), coef_names),
            centred   = setNames(coef(logit_interaction), coef_names)), 3)

# Likelihood-ratio test of the product term: a log-odds-scale question.
anova(logit, logit_interaction, test = "Chisq")

interaction_grid <- expand.grid(
  pratio = seq(price_range[[1]], price_range[[2]], length.out = 100),
  disp_coke = 0:1, disp_pepsi = 0)
interaction_plotdat <- as_tibble(predictions(logit_interaction,
                                             newdata = interaction_grid))
figure_plot <- ggplot(interaction_plotdat, aes(pratio, estimate,
                                colour = factor(disp_coke), fill = factor(disp_coke))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_rug(data = coke_in_range, aes(x = pratio), inherit.aes = FALSE,
           sides = "b", alpha = 0.1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(colour = "Coke display", fill = "Coke display",
       x = "Coke price / Pepsi price", y = "P(buy Coke)") + theme_minimal()
print(figure_plot)
ggsave("figures/interaction.png", figure_plot, width = 8, height = 4.5, dpi = 150)

# Price slopes at two prices under each display status.
slopes(logit_interaction, variables = "pratio", newdata = slope_grid)

# R8 decomposition by hand (centred parametrization: x enters as x - 1).
bi_all <- coef(logit_interaction)
bi <- c(alpha = bi_all[["(Intercept)"]],
        bP    = bi_all[["I(pratio - 1)"]],
        bC    = bi_all[["disp_coke"]],
        bPC   = bi_all[["I(pratio - 1):disp_coke"]])
decompose <- function(x) {
  e0 <- bi[["alpha"]] + bi[["bP"]] * (x - 1)
  e1 <- e0 + bi[["bC"]] + bi[["bPC"]] * (x - 1)
  c(pratio = x,
    product_term = bi[["bPC"]] * dlogis(e1),
    link_induced = bi[["bP"]] * (dlogis(e1) - dlogis(e0)),
    total        = (bi[["bP"]] + bi[["bPC"]]) * dlogis(e1) -
                   bi[["bP"]] * dlogis(e0))
}
rbind(decompose(0.9), decompose(1.3))

# Probability-scale contrasts with delta-method SEs, at both prices.
# Rows of slope_grid: (0.9, 0), (1.3, 0), (0.9, 1), (1.3, 1).
bind_rows(
  `pratio = 0.9` = as_tibble(slopes(logit_interaction, variables = "pratio",
                                    newdata = slope_grid, hypothesis = "b3 - b1 = 0")),
  `pratio = 1.3` = as_tibble(slopes(logit_interaction, variables = "pratio",
                                    newdata = slope_grid, hypothesis = "b4 - b2 = 0")),
  .id = "at") %>%
  select(at, estimate, std.error, conf.low, conf.high)

# ---- likelihood-geometry ----
X <- model.matrix(logit); y <- coke$coke

# Score at the MLE is zero: residuals are orthogonal to col(X).
crossprod(X, y - fitted(logit))
c(mean_fitted = mean(fitted(logit)), mean_y = mean(y))

# Newton-Raphson / IRLS from zero reproduces glm().
beta <- rep(0, ncol(X))
for (iter in 1:8) {
  eta_t <- drop(X %*% beta); p_t <- plogis(eta_t); w <- p_t * (1 - p_t)
  z_t   <- eta_t + (y - p_t) / w
  beta  <- drop(solve(crossprod(X, w * X), crossprod(X, w * z_t)))
}
w_hat <- fitted(logit) * (1 - fitted(logit))
cbind(IRLS = beta, glm = coef(logit),
      se_by_hand = sqrt(diag(solve(crossprod(X, w_hat * X)))),
      se_glm = sqrt(diag(vcov(logit))))

# ---- delta-method ----
V <- vcov(logit)
grad <- p0 * (1 - p0) * c(0, 1, 0, 0) +                     # f(eta0) e_k
        b[["pratio"]] * p0 * (1 - p0) * (1 - 2 * p0) * x0     # beta_k f'(eta0) x0
c(se_by_hand = sqrt(drop(t(grad) %*% V %*% grad)),
  se_package = slopes(logit, variables = "pratio", newdata = profile0)$std.error)

# ---- session-info ----
# Package versions matter: the `hypothesis` syntax in marginaleffects has changed over time.
sessionInfo()
