library(readxl)
library(tidyverse)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(progress)



set.seed(123)
# Prepare the data ----
# For now just focus on 1 study and 1 trait
# Use study_id = 9 and trait = development rate of egg
data <- read_excel("Data/vector-bite-extracted-maggie.xlsx") 

data <- data %>% 
  filter(study_id == 9 & trait == "development rate of egg") %>% 
  select(study_id, mean_temp, min_temp, max_temp, DTR, trait, trait_def, response, units) %>% 
  # If DTR = NA (constant temp), change to the 0, otherwise change it to the range of fluctuation
  mutate(DTR = ifelse(is.na(DTR), 0, max_temp - min_temp))

# Rename columns:  mean_temp -> temp; trait -> trait_name; response -> trait
colnames(data)[colnames(data) == "mean_temp"] <- "temp"
colnames(data)[colnames(data) == "trait"] <- "trait_name"
colnames(data)[colnames(data) == "response"] <- "trait"

# Change DTR to factor type
data$DTR <- as.factor(data$DTR)



## Plot the data
data %>%
  ggplot(aes(x = temp, y = trait, colour = DTR)) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2) +
  lims(x = c(5, 45), y = c(0,0.13)) +
  labs(y = "Development rate of egg (1/days)", x = paste("Temperature (ºC)")) +
  theme_bw()


data_const <- data %>%
  filter(DTR == 0)

# Fit TPC - Truncated normally-distributed Briere ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 15)
    cf.Tm ~ dunif(30, 45)
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
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 4 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
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
trait <- data_const$trait
N.obs <- length(trait)
temp <- data_const$temp

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
plot(trait ~ jitter(temp, 0.5), xlim = c(5, 45), ylim = c(0, 0.13), data = data_const, 
     ylab = "Development rate of egg (1/days)", xlab = expression(paste("Temperature (",degree,"C)")))

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
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(30, 60)
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
Temp.xs <- seq(0, 50, 0.1) # temperature gradient to calculate derived quantities over
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
trait <- data_const$trait
N.obs <- length(trait)
temp <- data_const$temp

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
plot(trait ~ jitter(temp, 0.5), xlim = c(10, 55), ylim = c(0, 0.13), data = data_const, 
     ylab = "Development rate of egg (1/days)", xlab = expression(paste("Temperature (",degree,"C)")))

lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


