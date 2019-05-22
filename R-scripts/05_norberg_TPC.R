library(tidyverse)
library(stringr)
library(rootSolve)
library(bbmle)
library(extrafont)
library(cowplot)


tdata <- read_csv("Data/ExtractedDataAllStudies_JB.csv")


tdata2 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime == 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, sep = "_"))

tdata_var <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime != 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, sep = "_")) %>% 
	filter(study_ID == "159")

unique_temps <- tdata2 %>% 
	group_by(curve.id) %>% 
	distinct(temp, .keep_all = TRUE) %>%
	tally() %>% 
	filter(n > 3)

tdata3 <- tdata2 %>% 
	filter(curve.id %in% unique_temps$curve.id)
tdata_var <- tdata_var %>% 
	filter(curve.id %in% unique_temps$curve.id) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = rate) 


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

tdata3 %>% 
	filter(study_ID == "159") %>% 
ggplot(aes(x = temp, y = rate)) + geom_point(shape = 1, size = 2, color = "grey") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	# geom_line(aes(x = temperature, y = predicted_rate), data = fits, color = "purple") +
	xlim(0, 40) +ylab("Rate") + xlab("Temperature (°C)") +
	facet_wrap( ~ curve.id, scales = "free_y") 
ggsave("figures/all_TPCs.png", width = 20, height = 20)


dat.full <- tdata3 %>% 
	filter(study_ID == "159") %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = rate) %>% 
	# filter(study_ID == 111) %>% 
	filter(temp_regime == 0) 
unique(dat.full$curve.id)

#### from Mridul's code, get the best fits for both of the TPCs
nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}

# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
z.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
w.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
a.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve

# Loop through all curve.id.list values to estimate parameters for all curves

for(i in 1:length(curve.id.list)){
	print(i)
	
	# Take a subset of the data corressponding to the ith curve.id.list value
	dat<-subset(dat.full,dat.full$curve.id==curve.id.list[i])
	
	# guess starting values for parameters 'z' and 'w'
	z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
	w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
	
	## This loop fits the model using a range of different starting guesses. We choose the best one using AIC. This helps find good solutions even if there are
	# convergence problems.
	# Starting estimates for parameters 'a' and 'b' use a plausible range but with broadly spaced estimates to speed up fitting. 
	avals<-seq(-0.2,1.2,0.1)		
	bvals<-seq(-0.2,0.3,0.05)
	mod.list<-list()
	AIC.list<-c()
	
	for(ia in 1:length(avals)){
		for(ib in 1:length(bvals)){
			a.guess<-avals[ia]
			b.guess<-bvals[ib]
			res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
								skip.hessian=TRUE,data=dat))
			if(class(res2)!="try-error"){
				mod.list<-append(mod.list,fit)
				AIC.list<-append(AIC.list,AIC(fit))
			}
		}
	}
	
	# Identify the best model from the list and save coefficients and R^2 values
	if(!is.null(AIC.list)){
		bestmodind<-which(AIC.list==min(AIC.list))
		if(length(bestmodind)>1){
			bestmodind<-sample(bestmodind,1)
		}
		bestmod<-mod.list[[bestmodind]]
		cfs<-coef(bestmod)
		expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
		rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
	}
	
	# If the quick fit yielded poor results (low R^2), try a more thorough search through parameter space
	if(rsqr<0.95){
		avals<-seq(-0.2,1.2,0.02)
		bvals<-seq(-0.2,0.3,0.02)
		mod.list<-list()
		AIC.list<-c()
		for(ia in 1:length(avals)){
			for(ib in 1:length(bvals)){
				a.guess<-avals[ia]
				b.guess<-bvals[ib]
				res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
									skip.hessian=TRUE,data=dat))
				if(class(res2)!="try-error"){
					mod.list<-append(mod.list,fit)
					AIC.list<-append(AIC.list,AIC(fit))
				}
			}
		}
		# Identify the best model from the list and save coefficients and R^2 values
		bestmodind<-which(AIC.list==min(AIC.list))
		if(length(bestmodind)>1){
			bestmodind<-sample(bestmodind,1)
		}
		
		bestmod<-mod.list[[bestmodind]]
		cfs<-coef(bestmod)
		expected <- nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
		rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
	}
	
	
	# Use the curve fit to find Topt and the estimated maximum growth rate (i.e. growth rate at Topt)
	grfunc<-function(x){
		-nbcurve(x,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
	}
	optinfo<-optim(c(x=cfs[[1]]),grfunc)
	opt<-optinfo$par[[1]]
	maxgrowth<- -optinfo$value
	
	#stash results		
	rsqr.list[i]<-rsqr
	z.list[i]<-cfs[[1]]
	w.list[i]<-cfs[[2]]
	a.list[i]<-cfs[[3]]
	b.list[i]<-cfs[[4]]
	s.list[i]<-cfs[[5]]
	topt.list[i]<-opt
	maxgrowth.list[i]<-maxgrowth
	n.list[i]<-length(dat$temperature)
}
# fits111 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant

fits159 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
# fits2 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
# write_csv(fits2, "data-processed/norberg-fits.csv")

fits4 <- bind_rows(fits3, fits168)

write_csv(fits4, "data-processed/norberg-fits-all.csv")
write_csv(fits159, "data-processed/norberg-fits-study159.csv")
write_csv(fits, "data-processed/norberg-fits-all-location.csv")
fits4 <- read_csv("data-processed/norberg-fits-all.csv")
### now make the plots!

nbcurve1<-function(x){
	res<-fits2$a.list[1]*exp(fits2$b.list[1]*x)*(1-((x-fits2$z.list[1])/(fits2$w.list[1]/2))^2)
	res
}


prediction_function <- function(df) {
	df <- df
	nbcurve1 <- function(x){
	res <- df$a.list[1]*exp(df$b.list[1]*x)*(1-((x-df$z.list[1])/(df$w.list[1]/2))^2)
	res
	}
	
	temps <- seq(0, 40, length = 150)
	
	predictions <- sapply(temps, nbcurve1)
	predicted_rate_nb <- data.frame(temperature = temps, predicted_rate = predictions)
	return(predicted_rate_nb)
	
}


fits_split <- fits159 %>% 
	filter(!is.na(topt.list)) %>% 
	split(.$curve.id.list)

df <- fits_split[[3]]

fits <- fits_split %>% 
	map_df(prediction_function, .id = "curve.id")

fits_above_zero <- fits %>% 
	filter(predicted_rate > -5, predicted_rate < 150)


nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}

temps <- seq(0, 40, length = 150)

predictions <- sapply(temps, nbcurve1)
predicted_rate_nb <- data.frame(temperature = temps, predicted_rate = predictions)



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p +
geom_point(aes(x = temperature, y = growth.rate), data = dat.full, shape = 1, size = 2, color = "grey") +
	geom_point(aes(x = temp, y = rate), data = tdata_var, shape = 1, size = 2, color = "cadetblue") +
	geom_point(aes(x = temp, y = perdicted_rate_var), data = all_159_varc, alpha = 0.5, size = 2, color = "purple") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	geom_line(aes(x = temperature, y = predicted_rate), data = fits_above_zero, color = "orange") +
	 xlim(9, 40) +ylab("Rate") + xlab("Temperature (°C)") +
	facet_wrap( ~ curve.id, scales = "free") 
ggsave("figures/norberg-tpc-fits-study159.png", width = 12, height = 8)

### let's make predictions for the variable regime

tdata_var %>% 
	View

fits159b <- fits159 %>% 
	rename(curve.id = curve.id.list)

all_159_var <- left_join(tdata_var, fits159b, by = "curve.id")

all_159_varc <- all_159_var %>% 
	mutate(predicted_rate_min = a.list*exp(b.list*min_temp)*(1-((min_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(predicted_rate_max = a.list*exp(b.list*max_temp)*(1-((max_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(perdicted_rate_var = (predicted_rate_min + predicted_rate_max)/2)
	

### let's look into study 111

tdata111 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	# mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	# mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp), (min_temp + max_temp)/2, mean_temp)) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>%
	# filter(temp_regime == 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = rate) %>% 
	mutate(curve.id = paste(study_ID, trait, sep = "_")) %>% 
	filter(study_ID == "111")

tdata111 %>% 
	# filter(temp_regime == 0) %>% 
	ggplot(aes(x = temperature, y = growth.rate, color = factor(temp_regime), shape = treatment_name_study)) + geom_point() +
	facet_wrap( ~ long)


temps18_111 <- read_csv("Data/study111-temps-18C-mean.csv")
temps25_111 <- read_csv("Data/study111-temps-25C-mean.csv")


all_111_temps <- bind_rows(temps18_111, temps25_111) %>% 
	mutate(study_ID = "111")

all_111_temps %>% 
	ggplot(aes(x = day, y = temperature, color = factor(mean_temperature))) + geom_line()


fits11125 <- fits111 %>%
	filter(curve.id.list == "111_fecundity_25.761") %>% 
	mutate(study_ID = "111")

fits11139 <- fits111 %>% 
	filter(curve.id.list == "111_fecundity_39.891") %>% 
	mutate(study_ID = "111")
	

a111_25 <- left_join(all_111_temps, fits11125)
a111_39 <- left_join(all_111_temps, fits11139)


all_111_preds_25 <- a111_25 %>% 
	mutate(predicted_rate = a.list*exp(b.list*temperature)*(1-((temperature-z.list)/(w.list/2))^2)) %>%
	mutate(predicted_rate_pos = ifelse(predicted_rate < 0, 0, predicted_rate)) %>% 
	group_by(mean_temperature, curve.id.list) %>% 
	summarise(mean_rate = mean(predicted_rate_pos)) %>% 
	rename(temperature = mean_temperature,
		   growth.rate = mean_rate) 

all_111_preds_39 <- a111_39 %>% 
	mutate(predicted_rate = a.list*exp(b.list*temperature)*(1-((temperature-z.list)/(w.list/2))^2)) %>%
	mutate(predicted_rate_pos = ifelse(predicted_rate < 0, 0, predicted_rate)) %>% 
	group_by(mean_temperature, curve.id.list) %>% 
	summarise(mean_rate = mean(predicted_rate_pos)) %>% 
	rename(temperature = mean_temperature,
		   growth.rate = mean_rate)

all_preds <- bind_rows(all_111_preds_25, all_111_preds_39) %>% 
	rename(curve.id = curve.id.list)


p + geom_point(aes(x = temperature, y = growth.rate), data = dat.full, shape = 1, size = 2, color = "grey") +
	# geom_point(aes(x = temperature, y = growth.rate), data = tdata_var, shape = 1, size = 2, color = "cadetblue") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	geom_line(aes(x = temperature, y = predicted_rate), data = fits_above_zero, color = "cadetblue") +
	xlim(0, 40) + ylab("Rate") + xlab("Temperature (°C)") + ylim(-3, 110) +
	geom_point(aes(x = mean_temp, y = rate, shape = treatment_name_study),
			   data = filter(tdata2, study_ID == "111", trait == "fecundity", temp_regime == 1))  +
	geom_point(aes(x = temperature, y = growth.rate), data = all_preds, color = "purple") +
	facet_wrap( ~ curve.id, scales = "free") 
ggsave("figures/study111-predictions-rate-summation.png", width = 12, height = 4)

#### study 168


study168_temps <- read_csv("Data/study168-temperature-fluctuations.csv")


study168_temps %>% 
	rename(fluctuation = flutation_type) %>% 
	ggplot(aes(x = hour, y = temperature, color = fluctuation, group = fluctuation)) + geom_line()








