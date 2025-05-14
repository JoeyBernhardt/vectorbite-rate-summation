## Rate summation project - study ID = 9

## Table of content:
##   
##    1. Prepare the data
##    2. Fitting TPC
##        A. Briere function (Truncated normally-distributed)
##        B. Quadratic function (Truncated normally-distributed)
##        C. Full Sharpe-Schoolfield model function
##    3. Performing rate summation
##        A. Extract TPC under constant temperature
##        B. Time series of the fluctuating regime
##        C. Rate summation calculations
##        D. Calculate the TPC posterior summary data
##    4. Plotting

###########
###### 0. Set-up workspace ----
###########

## Load Packages
library(readxl)
library(tidyverse)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(progress)
library(ggsci) # colour palettes
library(ggpubr)
library(progress)

## Load functions
source("R-scripts/00_Functions.R")


###########
###### 1. Prepare the data ----
###########

## Use study_id = 9 and trait = development rate of egg
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
plot <- data %>%
  ggplot(aes(x = temp, y = trait, colour = DTR)) +
  geom_point(size = 2) +
  #geom_hline(yintercept = 0, lty = 2) +
  scale_x_continuous(limits = c(0, 50)) + 
  scale_y_continuous(limits = c(-0.005, 0.13)) +
  scale_color_manual(name = "",
                     values = c("0" = "gray12", "10" = "firebrick1"),
                     breaks = c("0", "10"),
                     labels = c("Constant treatment", "Fluctuating treatment")) +
  labs(y = "Development rate of egg (1/days)", x = paste("Temperature (ºC)")) +
  theme_bw()

plot
# ggsave("figures/Lilian/data.png", plot, height = 6.21, width = 12)

# Create a new data frame for data from constant and fluctuating conditions
data_const <- data %>%
  filter(DTR == 0)
data_fluc <- data %>%
  filter(DTR != 0)

###########
###### 2. Fitting TPC ----
###########

##### A. Briere (Truncated normally-distributed) ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(30, 50)
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
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (110000 - 10000) / 100 ] * 5 = 5000
ni <- 110000 # number of iterations in each chain
nb <- 10000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 50, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

## Settings specific for fitting a normal/ truncated normal distribution ----

###  inits Function ----
inits<-function(){list(
  cf.q = runif(1, 0.01, 0.05),
  cf.T0 = runif(1, 5, 10),
  cf.Tm = runif(1, 35, 40),
  cf.sigma = rlnorm(1, log(1), 0.5)
  # cf.q = 0.01,
  # cf.T0 = 5,
  # cf.Tm = 40,
  # cf.sigma = rlnorm(1)
)
  }

### Parameters to Estimate ----
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

## Organize Data for JAGS ----
trait <- data_const$trait
N.obs <- length(trait)
temp <- data_const$temp

### define all the data for JAGS in a list object ----
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

## Run JAGS!! ----
set.seed(1234)
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

## Save the model
#save(model_briere, file = "R-scripts/R2jags-objects/study_9_mod13_briere.Rdata")

# Read the .rdata
load("R-scripts/R2jags-objects/study_9_mod13_briere.Rdata")

## Diagnostics ----
### Examine model output ----
model_briere$BUGSoutput$summary[1:5,]

### Run diagnostics ----
#mcmcplot(model_briere)

# Extract the DIC for future model comparisons
model_briere$BUGSoutput$DIC


## Plot the TPC ----
df.bri <- data.frame(model_briere$BUGSoutput$summary)
df.bri <- df.bri[-(1:5),]

# Add the corresponding temp to the dataframe
Temp.xs <- seq(0, 50, 0.1) # temperature gradient to calculate derived quantities over
df.bri$temp <- Temp.xs
df.bri <- df.bri %>% 
  relocate(temp)
df.bri$temp <- round(df.bri$temp, 1)

head(df.bri)

# Plot
bri.plot <- df.bri %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean, color = "Constant TPC"), linewidth = 1) +
  geom_point(data = data_const, aes(x = temp, y = trait, color = "Constant treatment"), , size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 50)) + 
  scale_y_continuous(limits = c(-0.005, 0.13)) +
  scale_color_manual(name = "",
                     values = c("Constant treatment" = "gray12", "Constant TPC" = "#000075"),
                     breaks = c("Constant treatment", "Constant TPC" )) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days^-1)",
    #title = "Briere model fit"
  ) +
  theme_bw()

# bri.plot <- df.bri %>% 
#   ggplot(aes(x = temp)) +
#   #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
#   geom_line(aes(y = mean), color = "black", linewidth = 1.5) +
#   # Customize the axes and labels
#   scale_x_continuous(limits = c(10, 45), breaks = seq(10, 46, 2)) + 
#   #scale_y_continuous(limits = c(-0.005, 0.15)) +
#   labs(
#     x = expression(paste("Temperature (", degree, "C)")),
#     y = "Performance"
#   ) +
#   theme_bw()

bri.plot
#ggsave("figures/Lilian/study_9_mod13_briere.png", bri.plot, height = 6.21, width = 12)


##### B. Quadratic function ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("quad_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 25)
    cf.Tm ~ dunif(30, 50)
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
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (110000 - 10000) / 100 ] * 5 = 5000
ni <- 110000 # number of iterations in each chain
nb <- 10000 # number of 'burn in' iterations to discard
nt <- 100 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 50, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

###  inits Function ----
inits<-function(){list(
  cf.q = runif(1, 0.01, 0.1),
  cf.T0 = runif(1, 0, 10),
  cf.Tm = runif(1, 35, 45),
  cf.sigma = rlnorm(1)
)}

### Parameters to Estimate ----
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

## Organize Data for JAGS ----
trait <- data_const$trait
N.obs <- length(trait)
temp <- data_const$temp

### define all the data for JAGS in a list object ----
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

## Run JAGS!! ----
#set.seed(1234)
model_quad <- jags(data = jag.data, 
                   inits = inits, 
                   parameters.to.save = parameters, 
                   model.file = "quad_T.txt",
                   n.thin = nt, 
                   n.chains = nc, 
                   n.burnin = nb, 
                   n.iter = ni, 
                   DIC = T, 
                   working.directory = getwd()
)

## Save the model
#save(model_quad, file = "R-scripts/R2jags-objects/study_9_mod10_quad.RData")

# Read the .rdata
load("R-scripts/R2jags-objects/study_9_mod9_quad.Rdata")

## Diagnostics ----
### Examine model output ----
head(model_quad$BUGSoutput$summary)

### Run diagnostics ----
#mcmcplot(model_quad)

# Extract the DIC for future model comparisons
model_quad$BUGSoutput$DIC


## Plot the TPC ----
df.quad <- data.frame(model_quad$BUGSoutput$summary)
df.quad <- df.quad[-(1:5),]

# Add the corresponding temp to the dataframe
df.quad$temp <- Temp.xs
df.quad <- df.quad %>% 
  relocate(temp)
df.quad$temp <- round(df.quad$temp, 1)

head(df.quad)

# Plot
quad.plot <- df.quad %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) +
  geom_line(aes(y = mean), size = 1) +
  geom_point(data = data_const, aes(x = temp, y = trait), color = "gray12", size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(10, 50)) + 
  #scale_y_continuous(limits = c(-0.005, 0.15)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days-1)",
    title = "Quadratic model fit"
  ) +
  theme_bw()

quad.plot
#ggsave("figures/Lilian/study_9_mod9_quad.png", quad.plot, height = 6.21, width = 9.29)


##### C. Full Sharpe-Schoolfield model function ----
## Model file ----
## write the model for JAGS and save it as a text file

sink("sharpe_full.txt")
cat("
    model{

    ## Priors
    r_tref ~ dunif(0, 10) # Rate at standardized temperature
    e ~ dunif(0, 1.5) # Activation energy (eV)
    el ~ dunif(1, 5) # low temperature de-activation energy (eV)
    tl ~ dunif(10, 30) # Temperature at which 1/2 enzyme suppressed due to low temperatures
    eh ~ dunif(1, 10) # High temperature de-activation energy (eV)
    th ~ dunif(30, 50) # Temperature at which 1/2 enzyme suppressed due to high temperatures
    sigma ~ dunif(0, 1000)
    tau <- 1 / (sigma * sigma)
    
    tref <- 20 # Standardisation temperature: rates are not inactivated by either high or low temperatures
    k <- 8.62e-05 # Boltzmann's constant

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- r_tref * exp(e/k * (1/(tref + 273.15) - 1/(temp[i] + 273.15))) / (1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp[i] + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(temp[i] + 273.15))))
    trait[i] ~ dnorm(trait.mu[i], tau)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- r_tref * exp(e/k * (1/(tref + 273.15) - 1/(Temp.xs[i] + 273.15))) / (1 + exp(-el/k * (1/(tl + 273.15) - 1/(Temp.xs[i] + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(Temp.xs[i] + 273.15))))
    }

    } # close model
    ",fill=T)
sink()

## MCMC Settings ----
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (1100000 - 100000) / 1000 ] * 5 = 5000
ni <- 1100000 # number of iterations in each chain
nb <- 100000 # number of 'burn in' iterations to discard
nt <- 1000 # thinning rate - jags saves every nt iterations in each chain
nc <- 5 # number of chains

## Derived Quantity Settings ----
Temp.xs <- seq(0, 50, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

###  inits Function ----
inits<-function(){list(
  r_tref = runif(1, 0, 10), # Rate at standardized temperature
  e = runif(1, 0, 0.5), # Activation energy (eV)
  el = runif(1, 1, 2.5), # low temperature de-activation energy (eV)
  tl = runif(1, 10, 20), # Temperature at which 1/2 enzyme suppressed due to low temperatures
  eh = runif(1, 1, 5), # High temperature de-activation energy (eV)
  th = runif(1, 30, 40), # Temperature at which 1/2 enzyme suppressed due to high temperatures
  sigma = rlnorm(1)
)}

### Parameters to Estimate ----
parameters <- c("r_tref", "e", "el", "tl", "eh", "th", "sigma", "z.trait.mu.pred")

## Organize Data for JAGS ----
trait <- data_const$trait
N.obs <- length(trait)
temp <- data_const$temp

### define all the data for JAGS in a list object ----
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

## Run JAGS!! ----
set.seed(1234)
model_sharpe <- jags(data = jag.data, 
                   inits = inits, 
                   parameters.to.save = parameters, 
                   model.file = "sharpe_full.txt",
                   n.thin = nt, 
                   n.chains = nc, 
                   n.burnin = nb, 
                   n.iter = ni, 
                   DIC = T, 
                   working.directory = getwd()
)

## Save the model
#save(model_sharpe, file = "R-scripts/R2jags-objects/study_9_mod13_sharpe.RData")

# Read the .rdata
load("R-scripts/R2jags-objects/study_9_mod12_sharpe.Rdata")

## Diagnostics ----
### Examine model output ----
model_sharpe$BUGSoutput$summary[1:8,]

### Run diagnostics ----
#mcmcplot(model_sharpe)

# Extract the DIC for future model comparisons
model_sharpe$BUGSoutput$DIC


## Plot the TPC ----
df.sharpe <- data.frame(model_sharpe$BUGSoutput$summary)
df.sharpe <- df.sharpe[-(1:8),]

# Add the corresponding temp to the dataframe
df.sharpe$temp <- Temp.xs
df.sharpe <- df.sharpe %>% 
  relocate(temp)
df.sharpe$temp <- round(df.sharpe$temp, 1)

head(df.sharpe)

# Plot
sharpe.plot <- df.sharpe %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "violet", alpha = 0.5) +
  geom_line(aes(y = mean), color = "darkorchid", size = 1) +
  geom_point(data = data_const, aes(x = temp, y = trait), color = "gray12", size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(5, 50)) + 
  #scale_y_continuous(limits = c(-0.005, 0.15)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days-1)",
    title = "Full Sharpe-Schoolfield model fit"
  ) +
  theme_bw()

sharpe.plot
#ggsave("figures/Lilian/study_9_mod11_sharpe.png", sharpe.plot, height = 6.21, width = 9.29)

all.plot <- ggarrange(bri.plot, quad.plot, sharpe.plot, nrow = 1, ncol = 3)
all.plot
#ggsave("figures/Lilian/study_9_allplot.png", all.plot, height = 6.21, width = 20)



###########
###### 3. Perform rate summation ----
###########

## For study 9, I chose the Briere function since it has a low DIC and the model
#fits better than Sharpe-Schoolfield

## Load the model
load("R-scripts/R2jags-objects/study_9_mod13_briere.Rdata")

## I will apply rate summation over temperature gradient 5ºC - 45ºC
## This is because DTR = 5ºC
temp_grad <- seq(0, 50, 0.1) # temperature gradient to apply rate summation. 


##### A. Extract TPC under constant temperature ----
## df: 501 columns (temperature gradient from 0 to 50ºC at 0.1ºC interval) x 5000 rows (MCMC iterations)
TPC_const <- data.frame(model_briere$BUGSoutput$sims.list$z.trait.mu.pred)



##### B. Time series of the fluctuating regime ----

## In study 9, individuals were placed in mean-5ºC for 24 hours then mean+5ºC for the next 24 hours and so on
## df = 501 cols (mean temperature gradient) x 2 rows
## Round to nearest 0.1ºC
timetemps_df <- t(data.frame(min = round(temp_grad - 5, 1), max = round(temp_grad + 5, 1)))
timetemps_df[timetemps_df < 0] <- 0
timetemps_df[timetemps_df > 50] <- 50
colnames(timetemps_df) <- temp_grad


##### C. Rate summation calculations ----
pred_RS <- RSCalcTempGrad(TPC_const, timetemps_df, temp_grad)

## Save the output of rate summation
#save(pred_RS, file = "data-processed/Lilian/study_9_predictions.Rdata")

## Load the output of rate summation
load("data-processed/Lilian/study_9_predictions.Rdata")


##### D. Calculate the TPC posterior summary data ----
pred_RS_summary <- calcPostQuants(pred_RS, "development rate", temp_grad)


###########
###### 4. Plotting ----
###########

RSplot <- pred_RS_summary %>% ggplot() +
  # Constant temp TPC
  geom_ribbon(data = df.bri, aes(x = temp, ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(data = df.bri, aes(x = temp, y = mean, color = "constant TPC"), linewidth = 1) +
  
  # Rate summation predicions
  geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI), fill = "gold1", alpha = 0.5) +
  geom_line(aes(x = temperature, y = mean, color = "RS predictions"), size = 1) +
  
  # Actual data from fluctuating temp
  geom_point(data = data_fluc, aes(x = temp, y = trait, color = "Fluctuating treatment"), size = 2) +
  
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(-0.005, 0.13)) +
  labs(
   x = expression(paste("Temperature (", degree, "C)")),
   y = "Development rate of egg (days-1)"
  ) +
  scale_color_manual(name = "",
                     values = c("Fluctuating treatment" = "firebrick1", "constant TPC" = "#000075", "RS predictions" = "#EFC000FF"),
                     breaks = c("Fluctuating treatment", "constant TPC", "RS predictions")) +
  theme_bw()

RSplot

## Save the plot
# ggsave("figures/Lilian/study_9_pred.png", RSplot, height = 6.21, width = 12)


