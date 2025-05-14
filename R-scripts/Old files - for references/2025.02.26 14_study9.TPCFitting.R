## Rate summation project - study ID = 9

## Table of content:
##   
##    1. Prepare the data
##    2. Fitting TPC
##        A. Briere function (Truncated normally-distributed)
##        B. Quadratic function (Truncated normally-distributed)
##        C. Full Sharpe-Schoolfield model function
##    3. Performing rate summation


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
data %>%
  ggplot(aes(x = temp, y = trait, colour = DTR)) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2) +
  lims(x = c(5, 45), y = c(0,0.13)) +
  labs(y = "Development rate of egg (1/days)", x = paste("Temperature (ºC)")) +
  theme_bw()

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
# plot(trait ~ jitter(temp, 0.5), xlim = c(5, 45), ylim = c(0, 0.13), data = data_const, 
#      ylab = "Development rate of egg (1/days)", xlab = expression(paste("Temperature (",degree,"C)")))
# 
# lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

### Use ggplot
df.bri <- data.frame(model_briere$BUGSoutput$summary)
df.bri <- df.bri[-(1:5),]

# Add the corresponding temp to the dataframe
df.bri$temp <- Temp.xs
df.bri <- df.bri %>% 
  relocate(temp)
df.bri$temp <- round(df.bri$temp, 1)

head(df.bri)

# Plot
bri.plot <- df.bri %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean), color = "#000075", linewidth = 1) +
  geom_point(data = data_const, aes(x = temp, y = trait), color = "gray12", size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(5, 50)) + 
  #scale_y_continuous(limits = c(-0.005, 0.15)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days^-1)",
    title = "Briere model fit"
  ) +
  theme_bw()

bri.plot
#ggsave("figures/Lilian/study_9_mod13_briere.png", bri.plot, height = 6.21, width = 9.29)


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

################################################################################
# Rate summation ----

# Define a function to calculate trait value under fluctuating temperature
rate_summation <- function(mean_temp, tpc){
  temp_max <- mean_temp + 5
  temp_min <- mean_temp - 5
  
  # Extract the trait values from the TPC 
  trait_max <- tpc %>%
    filter(temp == temp_max) %>%
    pull(mean)
  
  x2.5_max <- tpc %>%
    filter(temp == temp_max) %>%
    pull(X2.5.)
  
  x97.5_max <- tpc %>%
    filter(temp == temp_max) %>%
    pull(X97.5.)
  
  trait_min <- tpc %>% 
    filter(temp == temp_min) %>% 
    pull(mean)
  
  x2.5_min <- tpc %>%
    filter(temp == temp_min) %>%
    pull(X2.5.)
  
  x97.5_min <- tpc %>%
    filter(temp == temp_min) %>%
    pull(X97.5.)
  
  HPDI.upper <- mean((x97.5_min - trait_min), (x97.5_max - trait_max))
  HPDI.lower <- mean((trait_min - x2.5_min), (trait_max - x2.5_max))
  
  # Debugging outputs
  cat("Mean Temp:", mean_temp, "\n")
  cat("Temp Max:", temp_max, "Trait Max:", trait_max, "\n")
  cat("Temp Min:", temp_min, "Trait Min:", trait_min, "\n")
  
  # Calculate the average trait values for fluctuating temperature
  # If missing one value, then return NA
  fluct_trait <- c(mean(c(trait_max, trait_min), na.rm = F), HPDI.lower, HPDI.upper)
  cat("Trait value under fluctuating temp:", fluct_trait, "\n\n")
  return(fluct_trait)
}

# Predict the trait value under fluctuating temp
Temp_to_pred <- seq(5,40,0.5)
fluc_trait <- sapply(Temp_to_pred, rate_summation, tpc = df.bri)

# Put the results into a data frame

df <- data.frame()

for (i in 1:ncol(fluc_trait)){
  mn = fluc_trait[,i][1]
  Hl = fluc_trait[,i][2]
  Hu = fluc_trait[,i][3]  

  df <- rbind(df, c(mn,Hl,Hu))
}

df$temp <- Temp_to_pred
colnames(df) <- c("mean", "hpdi.low", "hpdi.up", "temp")

rate.sum.plot <- df.bri %>% 
  ggplot(aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  geom_line(aes(y = mean, color = "Constant"), size = 1) +
  geom_ribbon(data = df, aes(, ymin = mean - hpdi.low, ymax = mean + hpdi.up), fill = "gold1", alpha = 0.5) +
  geom_line(data = df, aes(y = mean, color = "Fluctuate"), size = 1) +
  geom_point(data = data, aes(x = temp, y = trait, color = DTR), size = 2) +
  # Customize the axes and labels
  #scale_x_continuous(limits = c(5, 45)) + 
  #scale_y_continuous(limits = c(-0.005, 0.15)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days-1)",
    #title = "Briere model fit"
  ) +
  scale_color_manual(name = "",
                     values = c("#0073C2FF", 
                                "#D55E00",
                                "#000075", 
                                "#EFC000FF"),
                     labels = c("DTR = 0", "DTR = 10", "Constant TPC", "Fluctuating TPC")) +
  #scale_color_manual(values = c("constant" = "#000075", "fluctuating" = "#EFC000FF")) +
  theme_bw()

rate.sum.plot
#ggsave("figures/Lilian/study_9_pred_notpc.png", rate.sum.plot, height = 6.21, width = 12)


rate.sum.plot.notpc  <- ggplot(data = df) +
  #geom_ribbon(data = df.bri, aes(ymin = X2.5., ymax = X97.5.), fill = "#4363d8", alpha = 0.5) +
  #geom_line(data = df.bri, aes(y = mean, color = "Constant"), size = 1) +
  geom_ribbon(data = df, aes(x = temp, ymin = mean - hpdi.low, ymax = mean + hpdi.up), fill = "gold1", alpha = 0.5) +
  geom_line(data = df, aes(x = temp, y = mean, color = "Estimated TPC"), size = 1) +
  geom_point(data = data_fluc, aes(x = temp, y = trait, color = DTR), size = 2) +
  # Customize the axes and labels
  scale_x_continuous(limits = c(0, 45)) + 
  #scale_y_continuous(limits = c(-0.005, 0.15)) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Development rate of egg (days-1)",
    #title = "Briere model fit"
  ) +
  scale_color_manual(name = "DTR",
                     values = c(#"#0073C2FF", 
                                "#D55E00",
                                #"#000075", 
                                "#EFC000FF"),
                     #labels = c("DTR = 0", "DTR = 10", "Constant T", "Fluctuating T")
                     ) +
  #scale_color_manual(values = c("constant" = "#000075", "fluctuating" = "#EFC000FF")) +
  theme_bw()

rate.sum.plot.notpc 
#ggsave("figures/Lilian/study_9_pred_notpc.png", rate.sum.plot.notpc, height = 6.21, width = 12)
