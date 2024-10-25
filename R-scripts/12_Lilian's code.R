# Load libraries
library(tidyverse)
library(nimble)
library(HDInterval)
library(MCMCvis)
library(coda) # makes diagnostic plots
library(IDPmisc) # makes nice colored pairs plots to look at joint posteriors
library(matrixStats)
library(truncnorm)
library(tidyverse)
library(janitor)
library(bayesTPC)

# For now just focus on 1 study and 1 trait ----
# Use Kingsolver growth rate data

king_const <- read_csv("Data/Kingolver-constant-25.csv") %>% 
  mutate(temperature = round(temperature, digits = 0)) %>% 
  mutate(growth_rate = exp(growth_rate)) %>% ## Since the raw data is ln(growth_rate)
  mutate(DTR = 0) %>% 
  mutate(temp_regime = 0) %>% ## Constant temp
  rename(Temp = temperature, Trait = growth_rate) ## Since bayesTPC expects data to be in a named list with the “Trait” as the response and “Temp” as the predictor

king_var <- read_csv("Data/Kingsolver-variable-growth.csv") %>% 
  mutate(mean_temp = round(mean_temp, digits = 0)) %>% 
  mutate(growth_rate_g = exp(growth_rate_g)) %>% ## Since the raw data is ln(growth_rate)
  mutate(temp_regime = 1) %>%  ## fluctuating temp
  rename(DTR = fluctuation, Trait = growth_rate_g, Temp = mean_temp) %>% 
  mutate(min_temp = Temp - DTR) %>% ## Calculate the min. temp
  mutate(max_temp = Temp + DTR) %>% ## Calculate the max. temp
  filter(DTR != 0) # Remove constant temp from this dataset since it doesn't match with king_const
  
## Combine the two datasets
king_comb <- bind_rows(king_const, king_var)

king_comb <- king_comb %>%
  # Change the min and max temp of the constant temp treatment into mean temp
  mutate(min_temp = ifelse(is.na(min_temp), Temp, min_temp)) %>%
  mutate(max_temp = ifelse(is.na(max_temp), Temp, max_temp))

# Plot the data
king_comb %>% 
  ggplot(aes(x = Temp, y = Trait, color = as.factor(DTR))) +
  geom_point()

# Fit TPC ----
set.seed(1234)

# Fit the Briere function
get_formula("briere")

# We only fit the TPC to the constant temp treatment
king_fit_briere <- b_TPC(data = king_const,
                     model = "briere",
                     niter = 11000, ## number of iteration
                     burn = 1000, ## number of burn in samples
                     samplerType = "AF_slice", ## slice sampler
                     priors = list(T_min = "dunif(-20, 20)",
                                   T_max = "dunif(30, 70)",
                                   sigma.sq = "dexp(1)") ## priors
                     )

summary(king_fit_briere)

# Look at the trace plot of all the parameters
par(mfrow = c(2,2))
traceplot(king_fit_briere, burn = 1000)

# Look at the ACF of the chains
s1 <- as.data.frame(king_fit_briere$samples[1000:10000,])
par(mfrow = c(2,2))
for(i in 1:4) {
  acf(s1[,i], lag.max = 50, main = "",
      ylab = paste("ACF: ", names(s1)[i], sep =""))
}

# Compare the prior and the posterior distributions
par(mfrow = c(1,1))
ppo_plot(king_fit_briere, burn = 1000)
