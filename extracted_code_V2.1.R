## ----------------------------------------------------------------------------------
#| warning: false

library(tidyverse)
library(broom)


## ----------------------------------------------------------------------------------
load(url("http://www.principlesofeconometrics.com/poe5/data/rdata/liquor5.rdata"))
head(liquor5)

dta <- liquor5 |>
  rename(i = hh, t = year, y = liquor, x = income)


## ----------------------------------------------------------------------------------
dta |>
  nest(data = -i) |>
  mutate(fit = map(data, ~ lm(y ~ x, data = .x)),
         tidied = map(fit, tidy)) |>
  unnest(tidied) |>
  filter(i == 1)


## ----------------------------------------------------------------------------------
dta |>
  nest(data = -i) |>
  mutate(fit = map(data, ~ lm(y ~ x, data = .x)),
         tidied = map(fit, tidy)) |>
  unnest(tidied) |>
  filter(term == "x") |>
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 40) +
  theme_bw()


## ----------------------------------------------------------------------------------
dtaSSW <- dta |>
  group_by(i) |>
  mutate(mx  = mean(x),
         my  = mean(y),
         Wxx = (x - mx)^2,
         Wxy = (x - mx) * (y - my),
         Wyy = (y - my)^2) |>
  ungroup()


## ----------------------------------------------------------------------------------
head(dtaSSW)


## ----------------------------------------------------------------------------------
hh1 <- filter(dtaSSW, i == 1)
b <- sum(hh1$Wxy) / sum(hh1$Wxx)
b


## ----------------------------------------------------------------------------------
a <- mean(hh1$y) - b * mean(hh1$x)
a


## ----------------------------------------------------------------------------------
tidy(lm(y ~ x, data = filter(dta, i == 1)))


## ----------------------------------------------------------------------------------
dtaSST <- dta |>
  mutate(gmx = mean(x),
         gmy = mean(y),
         Txx = (x - gmx)^2,
         Txy = (x - gmx) * (y - gmy),
         Tyy = (y - gmy)^2)


## ----------------------------------------------------------------------------------
head(dtaSST)


## ----------------------------------------------------------------------------------
b.pooled <- sum(dtaSST$Txy) / sum(dtaSST$Txx)
b.pooled


## ----------------------------------------------------------------------------------
a.pooled <- mean(dta$y) - b.pooled * mean(dta$x)
a.pooled


## ----------------------------------------------------------------------------------
tidy(lm(y ~ x, data = dta))


## ----------------------------------------------------------------------------------
#| warning: false

library(plm)

dtaPDF <- pdata.frame(dta, index = c("i", "t"))
fit.pooled <- plm(y ~ x, data = dtaPDF, model = "pooling")
tidy(fit.pooled)


## ----------------------------------------------------------------------------------
fit.dummy.variable <- lm(y ~ 0 + x + factor(i), data = dta)
head(tidy(fit.dummy.variable))


## ----------------------------------------------------------------------------------
fit.fixed <- plm(y ~ x, data = dtaPDF, model = "within")
tidy(fit.fixed)


## ----------------------------------------------------------------------------------
b.within <- sum(dtaSSW$Wxy) / sum(dtaSSW$Wxx)
b.within


## ----------------------------------------------------------------------------------
a.within <- mean(hh1$y) - b.within * mean(hh1$x)
a.within


## ----------------------------------------------------------------------------------
fixef(fit.fixed)[1]


## ----------------------------------------------------------------------------------
household.slopes <- dtaSSW |>
  group_by(i) |>
  summarise(Wxx_i = sum(Wxx),
            b_i   = sum(Wxy) / sum(Wxx)) |>
  mutate(weight = Wxx_i / sum(Wxx_i))

sum(household.slopes$weight * household.slopes$b_i)  # equals b.within
mean(household.slopes$b_i)                           # simple average differs


## ----------------------------------------------------------------------------------
fit.random <- plm(y ~ x, data = dtaPDF, model = "random")
summary(fit.random)


## ----------------------------------------------------------------------------------
fit.random$ercomp$sigma2          # the variance components

sigma2_e     <- fit.random$ercomp$sigma2[["idios"]]  # idiosyncratic error
sigma2_alpha <- fit.random$ercomp$sigma2[["id"]]     # household effect

TT <- pdim(dtaPDF)$nT$T           # number of time periods (TT, since T means TRUE)
TT

phi2 <- sigma2_e / (sigma2_e + TT * sigma2_alpha)   # weight on between variation
phi2

theta <- 1 - sqrt(phi2)           # quasi-demeaning parameter (POE5 eq. 15.34)
theta
fit.random$ercomp$theta           # same as plm's theta


## ----------------------------------------------------------------------------------
Wxy <- sum(dtaSSW$Wxy)
Wxx <- sum(dtaSSW$Wxx)

Txy <- sum(dtaSST$Txy)
Txx <- sum(dtaSST$Txx)

Bxy <- Txy - Wxy
Bxx <- Txx - Wxx

# b pooled OLS (phi2 = 1)
Txy / Txx

# b fixed effects or OLS with dummy variables (phi2 = 0)
Wxy / Wxx

# b between (OLS on household means)
b.between <- Bxy / Bxx
b.between

# b random effects GLS estimator
b.random <- (Wxy + phi2 * Bxy) / (Wxx + phi2 * Bxx)
b.random

# same as
coef(fit.random)[2]

# the intercept in the random effects model
a.random <- mean(dta$y) - b.random * mean(dta$x)
a.random

# same as
coef(fit.random)[1]


## ----------------------------------------------------------------------------------
dtaQD <- dtaSSW |>
  mutate(y_star = y - theta * my,
         x_star = x - theta * mx)

coef(lm(y_star ~ x_star, data = dtaQD))[2]

w <- Wxx / (Wxx + phi2 * Bxx)
w * b.within + (1 - w) * b.between


## ----------------------------------------------------------------------------------
#| fig-cap: "Magnitude of panel data estimators as a function of the between weight"
#| fig-height: 3.5

tibble(weight = seq(0, 1, by = 0.01)) |>
  mutate(slope = (Wxy + weight * Bxy) / (Wxx + weight * Bxx)) |>
  ggplot(aes(x = weight, y = slope)) +
  geom_line() +
  geom_hline(yintercept = b.between, linetype = "dashed") +
  annotate("point", x = c(0, phi2, 1),
           y = c(b.within, b.random, b.pooled), size = 2.5) +
  annotate("text", x = c(0, phi2, 1),
           y = c(b.within, b.random, b.pooled),
           label = c("Fixed effects", "Random effects", "Pooled OLS"),
           hjust = c(-0.15, -0.1, 1.1), vjust = c(0.5, 1.6, 1.6)) +
  annotate("text", x = 0.02, y = b.between, label = "Between",
           hjust = 0, vjust = -0.6) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = expression(phi^2 ~ "(weight on between variation)"),
       y = "Slope estimate") +
  theme_bw()

