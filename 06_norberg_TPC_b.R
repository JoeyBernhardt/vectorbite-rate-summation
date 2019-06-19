

library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)

nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}


### load Kingsolver data

tdata_king <- read_csv("Data/Kingolver-constant-25.csv")

sink("norberg.txt") 
cat("model{ 
	## Priors
	cf.z ~ dunif(0, 100)
	cf.w ~ dunif(5, 30)
	cf.a ~ dunif(-0.02, 1)
	cf.b ~ dunif(0, 1)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	## Likelihood 
	for(i in 1:N.obs){
	trait.mu[i] <- -1 * cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2) trait[i] ~ dnorm(trait.mu[i], cf.tau) } 
	## Derived Quantities and Predictions 
	for(i in 1:N.Temp.xs){ z.trait.mu.pred[i] <- -1 * cf.a * exp(cf.b * Temp.xs[i]) * (1-((Temp.xs[i]-cf.z)/(cf.w/2))^2) } 
	} # close model ",fill=T)
sink()


## parameters to estimate
parameters <- c("cf.z", "cf.w", "cf.a", "cf.b", "cf.sigma", "z.trait.mu.pred")


# Initial values for the parameters
inits<-function(){list(
	cf.z = 60,
	cf.w = 15,
	cf.a = 0,
	cf.b = 0.01,
	cf.sigma = rlnorm(1))}


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.2)
N.Temp.xs <-length(Temp.xs)

### Fitting the trait thermal response; Pull out data columns as vectors
data <- tdata_king %>%
	mutate(growth_rate = exp(growth_rate)) %>% 
	rename(trait = growth_rate,
		   T = temperature) 
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)


lf.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
			   model.file="norberg.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
			   n.iter=ni, DIC=T, working.directory=getwd())	

lf.fit.mcmc <- as.mcmc(lf.fit)	

b_params <- as.data.frame(lf.fit$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")


plot(lf.fit.mcmc[,c(1,3,4)])

plot(trait ~ T, xlim = c(0.2, 45), ylim = c(-1,2), data = data, ylab = "Growth rate", xlab = "Temperature")
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

predictions <- as.data.frame(lf.fit$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions) <- Temp.xs

predictions_sub <- predictions %>% 
	mutate(iteration = rownames(.)) %>% 
	select(iteration, everything()) %>% 
	sample_n(size = 100)


predictions_long <- predictions_sub %>% 
	gather(key = temperature, value = growth_rate, 2:227)

predictions_summary <- predictions %>% 
	mutate(iteration = rownames(.)) %>% 
	select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:227) %>% 
	group_by(temperature) %>% 
	summarise(q2.5=quantile(growth_rate, probs=0.025),
			  q97.5=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature))


predictions_long2 <- predictions_long %>% 
	mutate(temperature = as.numeric(temperature))

king_fits <- read_csv("data-processed/norberg-fits-kingolver.csv")

temperatures <- seq(0, 45, 0.2)

nbcurve_king<-function(temp,z,w,a,b){
	res<- king_fits$a.list[[1]]*exp(king_fits$b.list[[1]]*temp)*(1-((temp-king_fits$z.list[[1]])/(king_fits$w.list[[1]]/2))^2)
	res
}

king_preds <- data.frame(temperature = temperatures, growth_rate = sapply(temperatures, nbcurve_king))


nbcurve_kingb <- function(temp,z,w,a,b){
	res <- b_params$mean[b_params$term == "cf.a"]*exp(b_params$mean[b_params$term == "cf.b"]*temp)*(1-((temp- b_params$mean[b_params$term == "cf.z"])/(b_params$mean[b_params$term == "cf.w"]/2))^2)
	res
}

king_predsb <- data.frame(temperature = temperatures, growth_rate = sapply(temperatures, nbcurve_kingb))



	ggplot() + 
	geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T - 273.15, y = trait), data = data, color = "purple") + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = growth_rate), color = "pink", data = king_preds) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary, color = "orange") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary, color = "orange", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary, color = "orange", linetype = "dashed") +
	geom_line(aes(x = temperature - 273.15, y = mean), data = predictions_summary_briere, color = "blue") +
	geom_line(aes(x = temperature - 273.15, y = q2.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
	geom_line(aes(x = temperature - 273.15, y = q97.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
		
	ylab("Growth rate") + xlab("Temperature (°C)") +
	ylim(0, 3) + xlim(-10, 45)
ggsave("figures/kingsolver_tpcb_100.png", width = 10, height = 8)

mcmcplot(lf.fit)


## move on to Briere model

55+273
sink("briere.txt") 
cat("model{ 
	## Priors
	cf.q ~ dunif(0, 1)
	cf.T0 ~ dunif(253, 297)
	cf.Tm ~ dunif(298, 328)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	## Likelihood 
	for(i in 1:N.obs){
	trait.mu[i] <- cf.q * temp[i]*(temp[i] - cf.T0) * sqrt(abs(cf.Tm - temp[i])) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
	trait[i] ~ dnorm(trait.mu[i], cf.tau) } 
	## Derived Quantities and Predictions 
for(i in 1:N.Temp.xs){z.trait.mu.pred[i] <- cf.q * Temp.xs[i]*(Temp.xs[i] - cf.T0) * sqrt(abs(cf.Tm - Temp.xs[i])) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) } 
	} # close model ",fill=T)
sink()

parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")
273-10
# Initial values for the parameters
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 273+40,
	cf.T0 = 263,
	cf.sigma = rlnorm(1))}
# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
Temp.xs <- seq(-15, 50, 0.2) + 273.15
N.Temp.xs <-length(Temp.xs)

### Fitting the trait thermal response; Pull out data columns as vectors
data <- tdata_king %>%
	mutate(growth_rate = exp(growth_rate)) %>% 
	rename(trait = growth_rate,
		   T = temperature) %>% 
	mutate(T = T + 273.15)
	# this lets us reuse the same generic code: we only change this first line
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)


lf.fit_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
			   model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
			   n.iter=ni, DIC=T, working.directory=getwd())	

lf.fit.mcmc_briere <- as.mcmc(lf.fit_briere)


predictions_briere <- as.data.frame(lf.fit_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions_briere) <- Temp.xs

predictions_sub_briere <- predictions_briere %>% 
	mutate(iteration = rownames(.)) %>% 
	select(iteration, everything()) %>% 
	sample_n(size = 100)
View(predictions_sub_briere)
dim(predictions_sub_briere)

predictions_long_briere <- predictions_sub_briere %>% 
	gather(key = temperature, value = growth_rate, 2:327)

predictions_summary_briere <- predictions_briere %>% 
	mutate(iteration = rownames(.)) %>% 
	select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:327) %>% 
	group_by(temperature) %>% 
	summarise(q2.5=quantile(growth_rate, probs=0.025),
			  q97.5=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature))


predictions_long2_briere <- predictions_long_briere %>% 
	mutate(temperature = as.numeric(temperature))

mcmcplot(lf.fit_briere)

cf.T0 <- -9.07
cf.q <- 0.000526
cf.Tm <- 43.41



briere_curve <- function(temp, cf.T0 = -9.07, cf.q = 0.000526, cf.Tm = 43.41){
	res <- cf.q * temp*(temp - cf.T0) * sqrt(abs(cf.Tm - temp))*(temp <= cf.Tm)*(temp >= cf.T0)
	res
}

briere_params <- as.data.frame(lf.fit_briere$BUGSoutput$summary[1:3,]) %>% 
	rownames_to_column(var = "term")

temperatures <- seq(-20, 50)

briere_preds <- data.frame(T = temperatures, trait = sapply(temperatures, briere_curve))

briere_preds %>% 
	ggplot(aes(x = T, y = trait)) + geom_line() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = -9.07)

