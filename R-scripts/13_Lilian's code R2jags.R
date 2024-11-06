library(tidyverse)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(progress)

set.seed(123)
# Prepare the data ----
# For now just focus on 1 study and 1 trait
## Use Kingsolver growth rate data

king_const <- read_csv("Data/Kingolver-constant-25.csv") %>% 
  mutate(temperature = round(temperature, digits = 0)) %>% 
  mutate(growth_rate = exp(growth_rate)) %>% ## Since the raw data is ln(growth_rate)
  mutate(DTR = 0) %>% 
  mutate(temp_regime = 0) %>% ## Constant temp
  rename(temp = temperature, trait = growth_rate) ## Since bayesTPC expects data to be in a named list with the “Trait” as the response and “Temp” as the predictor

# king_var <- read_csv("Data/Kingsolver-variable-growth.csv") %>% 
#   mutate(mean_temp = round(mean_temp, digits = 0)) %>% 
#   mutate(growth_rate_g = exp(growth_rate_g)) %>% ## Since the raw data is ln(growth_rate)
#   mutate(temp_regime = 1) %>%  ## fluctuating temp
#   rename(DTR = fluctuation, trait = growth_rate_g, temp = mean_temp) %>% 
#   mutate(min_temp = temp - DTR) %>% ## Calculate the min. temp
#   mutate(max_temp = temp + DTR) %>% ## Calculate the max. temp
#   filter(DTR != 0) # Remove constant temp from this dataset since it doesn't match with king_const
# 
# ## Combine the two datasets
# king_comb <- bind_rows(king_const, king_var)
# 
# king_comb <- king_comb %>%
#   # Change the min and max temp of the constant temp treatment into mean temp
#   mutate(min_temp = ifelse(is.na(min_temp), temp, min_temp)) %>%
#   mutate(max_temp = ifelse(is.na(max_temp), temp, max_temp))

## Plot the data
# king_comb %>% 
#   ggplot(aes(x = temp, y = trait)) +
#   geom_point() +
#   labs(y = "Growth rate", x = paste("Temperature (ºC)")) +
#   lims(x = c(0, 60), y = c(1.1, 2)) +
#   theme_bw()

king_const %>% 
  ggplot(aes(x = temp, y = trait)) +
  geom_point() +
  labs(y = "Growth rate", x = paste("Temperature (ºC)")) +
  lims(x = c(0, 60), y = c(1.1, 2)) +
  theme_bw()


# Fit TPC - Truncated normally-distributed Briere ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 15)
    cf.Tm ~ dunif(40, 60)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()

## MCMC Settings ----
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 60, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

## Settings specific for fitting a normal/ truncated normal distribution ----

###  inits Function ----
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 40,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

### Parameters to Estimate ----
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

## Organize Data for JAGS ----
trait <- king_const$trait
N.obs <- length(trait)
temp <- king_const$temp

### define all the data for JAGS in a list object ----
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

## Run JAGS!! ----
model_briere <- jags(data = jag.data, 
                     inits=inits, 
                     parameters.to.save=parameters, 
                     model.file="briere_T.txt",
                     n.thin=nt, 
                     n.chains=nc, 
                     n.burnin=nb, 
                     n.iter=ni, 
                     DIC=T, 
                     working.directory=getwd()
                     )

## Diagnostics ----
### Examine model output ----
head(model_briere$BUGSoutput$summary)

### Run diagnostics ----
mcmcplot(model_briere)

# Extract the DIC for future model comparisons
model_briere$BUGSoutput$DIC


## Plot the TPC ----
plot(trait ~ jitter(temp, 0.5), xlim = c(0, 60), ylim = c(1.2, 3.2), data = king_const, 
     ylab = "growth rate", xlab = expression(paste("Temperature (",degree,"C)")))

lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# Fit TPC - Quadratic function ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("quad.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 15)
    cf.Tm ~ dunif(40, 60)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()

## MCMC Settings ----
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 60, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

## Settings specific for fitting a normal/ truncated normal distribution ----

###  inits Function ----
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 40,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

### Parameters to Estimate ----
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

## Organize Data for JAGS ----
trait <- king_const$trait
N.obs <- length(trait)
temp <- king_const$temp

### define all the data for JAGS in a list object ----
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

## Run JAGS!! ----
model_quad <- jags(data = jag.data, 
                   inits = inits, 
                   parameters.to.save = parameters, 
                   model.file = "quad.txt",
                   n.thin = nt, 
                   n.chains = nc, 
                   n.burnin = nb, 
                   n.iter = ni, 
                   DIC = T, 
                   working.directory = getwd()
)

## Diagnostics ----
### Examine model output ----
head(model_quad$BUGSoutput$summary)

### Run diagnostics ----
mcmcplot(model_quad)

# Extract the DIC for future model comparisons
model_quad$BUGSoutput$DIC


## Plot the TPC ----
plot(trait ~ jitter(temp, 0.5), xlim = c(0, 60), ylim = c(1.2,2.1), data = king_const, 
     ylab = "growth rate", xlab = expression(paste("Temperature (",degree,"C)")))

lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
