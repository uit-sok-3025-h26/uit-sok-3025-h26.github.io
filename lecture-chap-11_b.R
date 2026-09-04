# Lecture Chapter 11 B: Two-stage least squares in the Fulton Fish Market
#
options(digits = 4)

library(tidyverse)
require(ivreg) || {install.packages("ivreg") ; require(ivreg)}

#browseURL("http://www.principlesofeconometrics.com/poe5/data/def/fultonfish.def")

# fultonfish.def
# 
# obs:  111
# 
# date lprice quan lquan mon tue wed thu stormy mixed rainy cold totr diff change
# 
# date            date
# lprice          log(Price) of whiting( a common type of fish) per pound
# quan            Quantity of whiting sold, pounds
# lquan           log(Quantity)
# mon             Monday
# tue             Tuesday
# wed             Wednesday
# thu             Thursday
# stormy          High wind and waves
# mixed           Mixed wind and waves
# rainy           Rainy day on shore
# cold            Cold day on shore
# totr            Total received
# diff            Inventory change = totr-quan
# change          = 1 if diff large

load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/fultonfish.rdata"))
fish <- as_tibble(fultonfish)
rm(fultonfish)

fish

# lprice and lquan are already natural logarithms in the source dataset.

# -----------------------------------------------------------------------------
# 1. Is simultaneity economically plausible?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# If daily catch is predetermined, could demand be estimated by OLS? What role
# do dealer inventories play in deciding whether supply can respond to price?

fish %>%
  # 1. Use mutate first to create the explicit x-axis variable
  mutate(day = row_number()) %>% 
  # 2. Reshape only the relevant columns to long format
  pivot_longer(
    cols = c(lprice, lquan), 
    names_to = "variable", 
    values_to = "value"
  ) %>%
  # 3. Pipe directly into ggplot
  ggplot(aes(x = day, y = value, color = variable)) +
  geom_line() + 
  scale_color_manual(
    values = c("lprice" = "firebrick", "lquan" = "steelblue"),
    labels = c("Log price", "Log quantity")
  ) +
  labs(
    title = "Daily log price and log quantity",
    x = "Trading day",
    y = "Log value",
    color = NULL 
  ) +
  theme_classic() + 
  theme(
    legend.position = c(0.99, 0.99),          
    legend.justification = c("right", "top"), 
    legend.background = element_blank()       
  )

# Adjusted scale_y_continuous() allows us to plot two series with different scales on the same graph.
adjustment_factor <- 9

fish %>%
  # 1. Create the explicit x-axis variable
  mutate(day = row_number()) %>% 
  
  # 2. Pipe directly into ggplot (Notice we removed pivot_longer)
  ggplot(aes(x = day)) +
  
  # 3. Plot lprice normally (primary axis)
  geom_line(aes(y = lprice, color = "lprice")) +
  
  # 4. Plot lquan but subtract adjustment_factor so it drops down to hover around 0
  geom_line(aes(y = lquan - adjustment_factor, color = "lquan")) +
  
  # 5. Set up the two y-axes
  scale_y_continuous(
    name = "Log price", 
    # The secondary axis adds adjustment_factor back to the labels to show true quantity
    sec.axis = sec_axis(~ . + adjustment_factor, name = "Log quantity") 
  ) +
  
  # 6. Formatting colors and legend
  scale_color_manual(
    values = c("lprice" = "firebrick", "lquan" = "steelblue"),
    labels = c("Log price", "Log quantity")
  ) +
  labs(
    title = "Daily log price and log quantity",
    x = "Trading day",
    color = NULL 
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    
    # Optional: Match the axis text colors to the lines for easier reading
    axis.title.y = element_text(color = "firebrick"),
    axis.text.y = element_text(color = "firebrick"),
    axis.title.y.right = element_text(color = "steelblue"),
    axis.text.y.right = element_text(color = "steelblue")
  )

# Scatter plot
fish %>%
  ggplot(aes(x = lquan, y = lprice)) +
  geom_point(color = "steelblue", alpha = 0.6) + # alpha adds slight transparency
  labs(
    title = "Relationship between Log Quantity and Log Price",
    x = "Log quantity",
    y = "Log price"
  ) +
  theme_classic()

# Scatter plot with trading day gradient
fish %>%
  mutate(day = row_number()) %>%
  ggplot(aes(x = lquan, y = lprice, color = day)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c() + # A color scale that represents continuous numbers well
  labs(
    title = "Log Price vs. Log Quantity over Time",
    x = "Log quantity",
    y = "Log price",
    color = "Trading Day"
  ) +
  theme_classic()

# The graph documents variation; the institutional inventory argument supplies
# the reason to model price and quantity as jointly determined.

# -----------------------------------------------------------------------------
# 2. Does stormy weather shift equilibrium price?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Storms reduce supply. What sign should STORMY have in the price reduced form?

price_reduced_form <- lm(
  lprice ~ stormy + mon + tue + wed + thu,
  data = fish
)

summary(price_reduced_form)
coef(summary(price_reduced_form))

storm_coefficient <- coef(price_reduced_form)["stormy"]
storm_exact_percent <- 100 * (exp(storm_coefficient) - 1)

storm_coefficient
storm_exact_percent

stopifnot(abs(storm_coefficient - 0.346406) < 1e-5)
stopifnot(abs(coef(summary(price_reduced_form))["stormy", "t value"] - 4.638681) < 1e-5)

# -----------------------------------------------------------------------------
# 3. DELIBERATE FAILURE: count instruments but ignore their strength
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# The weekdays are excluded from supply. Does that fact alone make the supply
# equation empirically estimable?

price_restricted <- lm(lprice ~ stormy, data = fish)
weekday_test <- anova(price_restricted, price_reduced_form)
weekday_test

weekday_F <- weekday_test$F[2]
weekday_p <- weekday_test$`Pr(>F)`[2]

stopifnot(abs(weekday_F - 0.618762) < 1e-5)
stopifnot(abs(weekday_p - 0.650111) < 1e-5)

# Formal identification is not enough: the weekday demand shifters provide very
# little first-stage price variation, so supply is weakly identified in practice.

# -----------------------------------------------------------------------------
# 4. What does naive OLS say about demand?
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# OLS uses all observed price variation. Should its elasticity necessarily match
# the elasticity based only on storm-induced price variation?

ols_demand <- lm(
  lquan ~ lprice + mon + tue + wed + thu,
  data = fish
)
summary(ols_demand)

# -----------------------------------------------------------------------------
# 5. CORRECTION: estimate fish demand with dedicated 2SLS software
# -----------------------------------------------------------------------------

# Included exogenous variables appear on both sides of the vertical bar.
# Stormy is excluded from demand and instruments the endogenous lprice variable.
iv_demand <- ivreg::ivreg(
  lquan ~ lprice + mon + tue + wed + thu |
    stormy + mon + tue + wed + thu,
  data = fish
)

summary(iv_demand, diagnostics = TRUE)
coef(summary(iv_demand))

stopifnot(abs(coef(iv_demand)["lprice"] + 1.119417) < 1e-5)
stopifnot(abs(coef(summary(iv_demand))["lprice", "Std. Error"] - 0.428645) < 1e-5)

# A 1% storm-induced increase in price is associated with an estimated 1.12%
# decrease in quantity demanded, conditional on the instrument assumptions.

# -----------------------------------------------------------------------------
# 6. Stress-test the exclusion restriction
# -----------------------------------------------------------------------------

# PREDICTION STOP:
# Could poor weather on shore reduce restaurant attendance and therefore affect
# fish demand directly? Rainy and Cold address two observable versions of that
# concern, although no control set can prove instrument validity.

iv_demand_weather_controls <- ivreg::ivreg(
  lquan ~ lprice + mon + tue + wed + thu + rainy + cold |
    stormy + mon + tue + wed + thu + rainy + cold,
  data = fish
)

summary(iv_demand_weather_controls, diagnostics = TRUE)

comparison <- data.frame(
  model = c("OLS demand", "2SLS baseline", "2SLS + onshore weather controls"),
  price_elasticity = c(
    coef(ols_demand)["lprice"],
    coef(iv_demand)["lprice"],
    coef(iv_demand_weather_controls)["lprice"]
  ),
  standard_error = c(
    coef(summary(ols_demand))["lprice", "Std. Error"],
    coef(summary(iv_demand))["lprice", "Std. Error"],
    coef(summary(iv_demand_weather_controls))["lprice", "Std. Error"]
  )
)

print(comparison, row.names = FALSE)

# -----------------------------------------------------------------------------
# Optional demonstration: why the supply equation should not be trusted
# -----------------------------------------------------------------------------

iv_supply_weak <- ivreg::ivreg(
  lquan ~ lprice + stormy |
    mon + tue + wed + thu + stormy,
  data = fish
)

summary(iv_supply_weak, diagnostics = TRUE)

# This output is shown only as a diagnostic demonstration. Do not report the
# supply coefficient as a credible structural estimate: the excluded weekday
# instruments have a first-stage F statistic of only about 0.62.
