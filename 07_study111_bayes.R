

library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)


### study 111 soup to nuts

tdata <- read_csv("Data/ExtractedDataAllStudies_JB.csv")

tdata_111 <- tdata %>% 
	filter(study_ID %in% c(111)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime == 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, notes2, sep = "_")) %>% 
	mutate(study_trait = paste(study_ID, trait, sep = "_")) %>% 
	filter(grepl("25.76", long))



sink("norberg_lnorm.txt") 
cat("model{ 
	## Priors
	cf.z ~ dunif(0, 100)
	cf.w ~ dunif(5, 30)
	cf.a ~ dunif(-0.02, 20)
	cf.b ~ dunif(0, 1)
	# cf.sigma ~ dexp(0.001)
	cf.tau ~ dnorm(1000, 1/(250^2))
	cf.sigma <- sqrt(1/cf.tau)
	## Likelihood 
	for(i in 1:N.obs){
	trait.mu[i] <- cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2) 
	trait[i] ~ dlnorm(trait.mu[i], cf.tau) } 
	## Derived Quantities and Predictions 
	for(i in 1:N.Temp.xs){ z.trait.mu.pred[i] <- exp(cf.a * exp(cf.b * Temp.xs[i]) * (1-((Temp.xs[i]-cf.z)/(cf.w/2))^2)) } 
	} # close model ",fill=T)
sink()

sink("norberg_tnorm.txt") 
cat("model{ 
	## Priors
	cf.z ~ dunif(0, 100)
	cf.w ~ dunif(5, 30)
	cf.a ~ dunif(-0.02, 20)
	cf.b ~ dunif(0, 1)
	cf.tau ~ dnorm(1000, 1/(250^2))
	cf.sigma <- sqrt(1/cf.tau)
	# cf.tau <- 1/(cf.sigma^2)
	# cf.sigma ~ dexp(0.0001)
	
	## Likelihood 
	for(i in 1:N.obs){
	trait.temp[i] <- cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2)
	trait.mu[i] <- trait.temp[i]*((cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2)) >0) 
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) } 
	## Derived Quantities and Predictions

	for(i in 1:N.Temp.xs){ 
	temporary[i] <- cf.a * exp(cf.b * Temp.xs[i]) * (1-((Temp.xs[i]-cf.z)/(cf.w/2))^2)	
	z.trait.mu.pred[i] <- cf.a * exp(cf.b * Temp.xs[i]) * (1-((Temp.xs[i]-cf.z)/(cf.w/2))^2)*(temporary[i] > 0) } 
	} # close model ",fill=T)
sink()


## parameters to estimate
parameters <- c("cf.z", "cf.w", "cf.a", "cf.b", "cf.tau", "z.trait.mu.pred")


# Initial values for the parameters
inits<-function(){list(
	cf.z = 25,
	cf.w = 22,
	cf.a = 12,
	cf.b = 0.01,
	cf.tau = 100
	# cf.sigma = .1
	)}


# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.2)
N.Temp.xs <-length(Temp.xs)

### Fitting the trait thermal response; Pull out data columns as vectors
data <- tdata_111 %>%
	select(rate, mean_temp) %>% 
	rename(trait = rate,
		   T = mean_temp)
# data <- bind_rows(data, data)


trait <- data$trait
N.obs <- length(trait)
temp <- data$T

### plot it

data %>% 
	ggplot(aes(x = T, y = trait)) + geom_point()


# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)


lf.fit_exp <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
			   model.file="norberg_tnorm.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
			   n.iter=ni, DIC=T, working.directory=getwd())	

lf.fit.mcmc <- as.mcmc(lf.fit_exp)	

b_params <- as.data.frame(lf.fit$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")


plot(lf.fit.mcmc[,c(1,3,4)])

plot(trait ~ T, xlim = c(0.2, 45), ylim = c(0,90), data = data, ylab = "Growth rate", xlab = "Temperature")
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.fit$BUGSoutput$summary[7:(7 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



predictions <- as.data.frame(lf.fit_exp$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
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


ggplot() + 
	geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T, y = trait), data = data, color = "purple") + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary, color = "orange") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary, color = "orange", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary, color = "orange", linetype = "dashed") +
	# geom_line(aes(x = temperature - 273.15, y = mean), data = predictions_summary_briere, color = "blue") +
	# geom_line(aes(x = temperature - 273.15, y = q2.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
	# geom_line(aes(x = temperature - 273.15, y = q97.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
	
	ylab("Growth rate") + xlab("Temperature (°C)") +
	ylim(0, 100) + xlim(0, 45) 
ggsave("figures/EN-study111-data-lognormal-exp-truncated-normal-error.png", width = 6, height = 4)

mcmcplot(lf.fit_exp, filename = "study-111-data-tnorm")


print(1/sqrt((0.05)))
dev.off()

### are the parameters correlated?