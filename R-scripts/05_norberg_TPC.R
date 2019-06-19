library(tidyverse)
library(stringr)
library(rootSolve)
library(bbmle)
library(extrafont)
library(cowplot)


tdata <- read_csv("Data/ExtractedDataAllStudies_JB.csv")
tdata_king <- read_csv("Data/Kingolver-constant-25.csv")


tdata2 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime == 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, notes2, sep = "_")) %>% 
	mutate(study_trait = paste(study_ID, trait, sep = "_"))

tdata_var119 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime != 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, notes2, sep = "_")) %>% 
	filter(study_ID == "119", trait == "weight gain") 

unique_temps <- tdata2 %>% 
	group_by(curve.id) %>% 
	distinct(temp, .keep_all = TRUE) %>%
	tally() %>% 
	filter(n > 3)

tdata3 <- tdata2 %>% 
	filter(curve.id %in% unique_temps$curve.id)
tdata_var_119 <- tdata_var119 %>% 
	filter(curve.id %in% unique_temps$curve.id) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = rate) 


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

tdata3 %>% 
	filter(study_ID == "119") %>% 
ggplot(aes(x = temp, y = rate)) + geom_point(shape = 1, size = 2, color = "grey") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	# geom_line(aes(x = temperature, y = predicted_rate), data = fits, color = "purple") +
	# xlim(0, 40) +
	ylab("Rate") + xlab("Temperature (°C)") +
	facet_wrap( ~ curve.id, scales = "free_y") 
ggsave("figures/all_TPCs.png", width = 20, height = 20)


dat.full <- tdata_king %>% 
	mutate(growth_rate = exp(growth_rate)) %>% 
	rename(growth.rate = growth_rate) %>% 
	mutate(curve.id = "1")

#### norberg curve
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
fits_king <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
fits119 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant

fits150 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant

fits168 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant

fits159 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
# fits2 <- data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
# write_csv(fits2, "data-processed/norberg-fits.csv")

fits4 <- bind_rows(fits3, fits168)

write_csv(fits4, "data-processed/norberg-fits-all.csv")
write_csv(fits_king, "data-processed/norberg-fits-kingolver.csv")
write_csv(fits159, "data-processed/norberg-fits-study159.csv")
write_csv(fits119, "data-processed/norberg-fits-study119.csv")
write_csv(fits168, "data-processed/norberg-fits-study168.csv")
write_csv(fits150, "data-processed/norberg-fits-study150.csv")
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


fits_split <- fits119 %>% 
	filter(!is.na(topt.list)) %>% 
	split(.$curve.id.list)

df <- fits_split[[3]]

fits <- fits_split %>% 
	map_df(prediction_function, .id = "curve.id")

fits_above_zero <- fits %>% 
	filter(predicted_rate > 0, predicted_rate < 150)


nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

tdata_var168 <- left_join(tdata_var, study168)
write_csv(tdata_var168, "data-processed/tdatavar168.csv")

p +
geom_point(aes(x = temperature, y = growth.rate), data = dat.full, shape = 1, size = 2, color = "grey") +
	geom_point(aes(x = temp, y = rate), data = tdata_var119, shape = 1, size = 2, color = "cadetblue") +
	geom_point(aes(x = temp, y = perdicted_rate_var), data = all_119_varc, alpha = 0.5, size = 2, color = "purple") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_point(aes(x = mean_temp, y = mean_performance), data = study150t2, alpha = 0.5, size = 2, color = "purple") +
	# geom_point(aes(x = min_temp, y = mean_performance), data = study150t2, alpha = 0.5, size = 2, color = "purple") +
	# geom_point(aes(x = max_temp, y = mean_performance), data = study150t2, alpha = 0.5, size = 2, color = "purple") +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	geom_line(aes(x = temperature, y = predicted_rate), data = fits_above_zero, color = "orange") +
	 xlim(9, 40) +ylab("Rate") + xlab("Temperature (°C)") +
	facet_wrap( ~ curve.id, scales = "free") 
ggsave("figures/norberg-tpc-fits-study119.png", width = 12, height = 8)
ggsave("figures/norberg-tpc-fits-study150.png", width = 12, height = 8)
ggsave("figures/norberg-tpc-fits-study159.png", width = 12, height = 8)
ggsave("figures/norberg-tpc-fits-study168.png", width = 12, height = 8)

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


### study 119
fits119b <- fits119 %>% 
	rename(curve.id = curve.id.list)

all_119_var <- left_join(tdata_var119, fits119b, by = "curve.id")

all_119_varc <- all_119_var %>% 
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


fits <- read_csv("data-processed/norberg-fits-all-location.csv") %>% 
	filter(curve.id.list %in% c("111_fecundity_25.761", "111_fecundity_39.891")) %>% 
	mutate(study_ID = "111")

fits11125 <- fits %>%
	filter(curve.id.list == "111_fecundity_25.761") %>% 
	mutate(study_ID = "111")

fits11139 <- fits %>% 
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

all_preds_111 <- bind_rows(all_111_preds_25, all_111_preds_39) %>% 
	rename(curve.id = curve.id.list)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
	# geom_point(aes(x = temperature, y = growth.rate), data = dat.full, shape = 1, size = 2, color = "grey") +
	# geom_point(aes(x = temperature, y = growth.rate), data = tdata_var, shape = 1, size = 2, color = "cadetblue") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	# geom_line(aes(x = temperature, y = predicted_rate), data = fits_above_zero, color = "cadetblue") +
	xlim(0, 40) + ylab("Rate") + xlab("Temperature (°C)") + ylim(-3, 110) +
	geom_point(aes(x = mean_temp, y = rate, shape = treatment_name_study),
			   data = filter(tdata2, study_ID == "111", trait == "fecundity", temp_regime == 1))  +
	# geom_point(aes(x = temperature, y = growth.rate), data = all_preds, color = "purple") +
	facet_wrap( ~ curve.id, scales = "free") 
ggsave("figures/study111-predictions-rate-summation.png", width = 12, height = 4)

#### study 168


study168_temps <- read_csv("Data/study168-temperature-fluctuations.csv")


study168_temps %>% 
	rename(fluctuation = flutation_type) %>% 
	mutate(realized_temp = temperature + min_temp) %>% 
	mutate(unique_regime = paste(fluctuation, min_temp, sep = "_")) %>%  
	ggplot(aes(x = hour, y = realized_temp, color = unique_regime, group = unique_regime)) + geom_line()
	

study168_tempsb <- study168_temps %>% 
	rename(fluctuation = flutation_type) %>% 
	mutate(realized_temp = temperature + min_temp) %>% 
	mutate(unique_regime = paste(fluctuation, min_temp, sep = "_")) %>% 
	mutate(curve.id.list = fits168$curve.id.list[[1]])

study168_tempsb %>% 
	ggplot(aes(x = hour, y = temperature, color = fluctuation, group = fluctuation)) + geom_line()




study168 <- left_join(study168_tempsb, fits168) %>% 
	group_by(unique_regime, DTR, min_temp) %>% 
	sample_n(size = 23, replace = FALSE) %>% 
	mutate(predicted_performance = a.list*exp(b.list*realized_temp)*(1-((realized_temp-z.list)/(w.list/2))^2)) %>% 
	group_by(unique_regime, DTR, min_temp) %>% 
	summarise(mean_performance = mean(predicted_performance),
			  mean_temperature = mean(realized_temp))

write_csv(study168, "data-processed/study168-predicted-variable.csv")


study168 %>% 
	group_by(unique_regime) %>% 
	distinct(realized_temp) %>% 
	tally()

study168 <- read_csv("data-processed/study168-predicted-variable.csv")

study150_temps <- read_csv("Data/study150-temperature-fluctuations.csv") %>% 
	mutate(curve.id.list = "150_Duration of incubation_sinensis_120.1551_NA")

study150t2 <- left_join(study150_temps, fits150) %>% 
	mutate(predicted_performance = a.list*exp(b.list*temperature)*(1-((temperature-z.list)/(w.list/2))^2)) %>% 
	summarise(mean_performance = mean(predicted_performance),
			  mean_temp = mean(temperature),
			  max_temp = max(temperature),
			  min_temp = min(temperature))
	
	

### plot all together

fits159 <- read_csv("data-processed/norberg-fits-study159.csv") %>% 
	mutate(study_ID = "159")
fits119 <- read_csv("data-processed/norberg-fits-study119.csv") %>% 
	mutate(study_ID = "119")
fits168 <- read_csv("data-processed/norberg-fits-study168.csv") %>% 
	mutate(study_ID = "168")
fits150 <- read_csv("data-processed/norberg-fits-study150.csv") %>% 
	mutate(study_ID = "150")
fits <- read_csv("data-processed/norberg-fits-all-location.csv") %>% 
	filter(curve.id.list %in% c("111_fecundity_25.761", "111_fecundity_39.891")) %>% 
	mutate(study_ID = "111")

fits_all <- bind_rows(fits, fits159, fits119, fits168, fits150) %>% 
	mutate(curve.id.list = case_when(curve.id.list == "111_fecundity_25.761" ~ "111_fecundity_Drosophila_25.761_0_Begin Miami flies",
								curve.id.list == "111_fecundity_39.891" ~ "111_fecundity_Drosophila_39.891_0_Begin New Jersey flies",
								TRUE ~ curve.id.list))



data_sel <- tdata2 %>% 
	# dplyr::filter(study_trait %in% c(study_trait)) %>% 
	dplyr::filter(study_trait %in% c("111_fecundity", "159_thorax length_melanogaster",
									  "159_thorax length",
									  "159_winglength",
									  "119_weight gain",
									  "168_fecundity/body weight",
									  "150_Duration of incubation")) %>% 
	mutate(curve.id.orig = curve.id) %>% 
	mutate(curve.id == ifelse(curve.id == "159_thorax length_melanogaster_-4.77913_F_see notes file for calculated equivalent developmental temperatures for fluctuating regimes, expt design details and variance specifics.",
							  "159_thorax length_melanogaster_-4.77913_F_NA", curve.id)) %>% 
	dplyr::mutate(curve.id == case_when(curve.id == "111_fecundity_Drosophila_25.761_0_Begin Miami flies" ~ 
								 	"111_fecundity_Drosophila_25.761",
								 curve.id == "111_fecundity_Drosophila_39.891_0_Begin New Jersey flies" ~ 
								 	"111_fecundity_Drosophila_39.891",
								 TRUE ~ curve.id))

data_sel2 <- data_sel %>% 
	filter(curve.id %in% c(pred_c$curve.id.new)) %>% 
	mutate(curve.id.new = curve.id)

unique(data_sel2$curve.id)

# write_csv(data_sel, "data-processed/tdata-sel.csv")
data_sel <- read_csv("data-processed/tdata-sel.csv")
pred_c <- read_csv("data-processed/pred_c.csv") 


length(unique(data_sel2$curve.id.new))
length(unique(pred_c$curve.id.new))

setdiff(unique(data_sel2$curve.id.new), unique(predictions2$curve.id.new))
setdiff(unique(data_sel2$curve.id.new), unique(pred_c$curve.id.new))

unique(predictions$curve.id)

predictions2 <- left_join(predictions, pred_c) %>% 
	filter(is.na(curve.id.new)) %>% View
	mutate(curve.id.new = ifelse(is.na(curve.id.new), "119_weight gain_niloticus_NA", curve.id.new))


write_csv(pred_c, "data-processed/pred_c.csv")

unique(fits_all$curve.id.list)

fsplit <- fits_all %>% 
	split(.$curve.id.list)



predictions <- fsplit %>% 
	map_df(prediction_function, .id = "curve.id")


predictions2 <- left_join(predictions, pred_c) %>% 
	# filter(is.na(curve.id.new)) %>% 
	mutate(curve.id.new = ifelse(is.na(curve.id.new), "119_weight gain_niloticus_NA", curve.id.new))

write_csv(predictions2, "data-processed/predictions2.csv")

predictions2 %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = mean_temp, y = rate), data = data_sel2) +
	facet_wrap( ~ curve.id.new, scales = "free")
	
curves_pred <- unique(predictions$curve.id)
data_sel_curves <- unique(data_sel$curve.id)
# 
# pred_c <- predictions %>% 
# 	distinct(curve.id)
# 
# library(fuzzyjoin)
# joined <- fuzzyjoin::stringdist_inner_join(pred_c, data_sel, by = "curve.id", max_dist = 5) %>% 
# 	select(curve.id.x, curve.id.y, everything()) 



# now bring in the predicted rates ----------------------------------------


#### study 168
study168 <- read_csv("data-processed/study168-predicted-variable.csv")


#### study 150
study150_temps <- read_csv("Data/study150-temperature-fluctuations.csv") %>% 
	mutate(curve.id.list = "150_Duration of incubation_sinensis_120.1551_NA")

study150t2 <- left_join(study150_temps, fits150) %>% 
	mutate(predicted_performance = a.list*exp(b.list*temperature)*(1-((temperature-z.list)/(w.list/2))^2)) %>% 
	summarise(mean_performance = mean(predicted_performance),
			  mean_temp = mean(temperature),
			  max_temp = max(temperature),
			  min_temp = min(temperature))

#### study 111

fits <- read_csv("data-processed/norberg-fits-all-location.csv") %>% 
	filter(curve.id.list %in% c("111_fecundity_25.761", "111_fecundity_39.891")) %>% 
	mutate(study_ID = "111")

fits11125 <- fits %>%
	filter(curve.id.list == "111_fecundity_25.761") %>% 
	mutate(study_ID = "111")

fits11139 <- fits %>% 
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

all_preds_111 <- bind_rows(all_111_preds_25, all_111_preds_39) %>% 
	rename(curve.id = curve.id.list)



#### study 119

fits119b <- read_csv("data-processed/norberg-fits-study119.csv") %>% 
	mutate(study_ID = "119") %>% 
	rename(curve.id = curve.id.list)

fits119b <- fits119 %>% 
	rename(curve.id = curve.id.list)



tdata_var119 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime != 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, notes2, sep = "_")) %>% 
	filter(study_ID == "119", trait == "weight gain") 

all_119_var <- left_join(tdata_var119, fits119b, by = "curve.id")

all_119_varc <- all_119_var %>% 
	mutate(predicted_rate_min = a.list*exp(b.list*min_temp)*(1-((min_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(predicted_rate_max = a.list*exp(b.list*max_temp)*(1-((max_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(predicted_rate_var = (predicted_rate_min + predicted_rate_max)/2)


### study 159

fits159b <- fits159 %>% 
	rename(curve.id = curve.id.list)


tdata_var159 <- tdata %>% 
	# filter(study_ID %in% c(168, 165, 159, 150, 124, 119, 111, 91, 67, 61, 58)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime != 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, notes2, sep = "_")) %>% 
	filter(study_ID == "159") %>% 
	mutate(study_ID = as.character(study_ID)) %>% 
	mutate(curve.id = str_replace(curve.id, "_NA", ""))

all_159_var <- left_join(tdata_var159, fits159b)

all_159_varc <- all_159_var %>% 
	mutate(predicted_rate_min = a.list*exp(b.list*min_temp)*(1-((min_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(predicted_rate_max = a.list*exp(b.list*max_temp)*(1-((max_temp-z.list)/(w.list/2))^2)) %>% 
	mutate(predicted_rate_var = (predicted_rate_min + predicted_rate_max)/2)

unique(all_159_varc$curve.id)

predictions2 %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = mean_temp, y = rate), data = data_sel2) +
	facet_wrap( ~ curve.id.new, scales = "free")

unique(predictions2$curve.id.new)



# now gather the variable predictions -------------------------------------


all159 <- all_159_varc %>% 
	select(curve.id, temp, predicted_rate_var) %>% 
	rename(temperature = temp) 
all119 <- all_119_varc %>% 
	select(curve.id, temp, predicted_rate_var) %>% 
	rename(temperature = temp)
	
all150 <- study150t2 %>% 
	select(mean_performance, mean_temp) %>% 
	rename(temperature = mean_temp) %>% 
	mutate(curve.id = "150_Duration of incubation_sinensis_120.1551_NA") %>% 
	rename(predicted_rate_var = mean_performance)

all111 <- all_preds_111 %>% 
	rename(predicted_rate_var = growth.rate)
	

all168 <- study168 %>% 
	rename(predicted_rate_var = mean_performance) %>% 
	rename(temperature = mean_temperature) %>% 
	mutate(curve.id = fits168$curve.id.list[[1]]) %>% 
	select(curve.id, temperature, predicted_rate_var)
	

all_var <- bind_rows(all168, all111, all159, all150, all119) %>% 
	rename(curve.id.new = curve.id) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "111_fecundity_25.761", "111_fecundity_Drosophila_25.761")) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "111_fecundity_39.891", "111_fecundity_Drosophila_39.891")) %>%
	mutate(curve.id.new = str_replace(curve.id.new, "150_Duration of incubation_sinensis_120.1551_NA", "150_Duration of incubation_sinensis_120.1551")) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "119_weight gain_niloticus_NA_NA_NA", "119_weight gain_niloticus_NA")) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "168_fecundity/body weight_n. nevadensis_116.4236111_NA", "168_fecundity/body weight_n. nevadensis_116.4236111")) %>% 
mutate(curve.id.new = str_replace(curve.id.new, "119_weight gain_stratiotes_NA_NA_NA", "119_weight gain_stratiotes_NA")) 

	

predictions2 %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = mean_temp, y = rate), data = data_sel2) +
	geom_point(aes(x = temperature, y = predicted_rate_var), data = all_var, color = "blue") +
	facet_wrap( ~ curve.id.new, scales = "free")

unique(all_var$curve.id.new)
unique(predictions2$curve.id.new)

setdiff(unique(all_var$curve.id.new), unique(predictions2$curve.id.new))




# now let’s gather the variable observed ----------------------------------

tdata_var <- tdata %>% 
	# filter(study_ID %in% c(168, 159, 150, 119, 111)) %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(temp_regime != 0) %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	mutate(curve.id = paste(study_ID, trait, species, long, notes, sep = "_")) %>% 
	filter(trait %in% c("weight gain", "Duration of incubation", "fecundity", "fecundity/body weight", "thorax length",
						"winglength")) %>% 
	filter(curve.id != "168_fecundity_n. nevadensis_116.4236111_NA") %>% 
	mutate(curve.id = str_replace(curve.id, "_0", "")) %>% 
	mutate(curve.id = str_replace(curve.id, "_NA_NA", "_NA")) %>% 
	mutate(curve.id = str_replace(curve.id, "nevadensis_116.4236111_NA", "nevadensis_116.4236111")) %>% 
	mutate(curve.id = str_replace(curve.id, "sinensis_120.1551_NA", "sinensis_120.1551")) %>% 
	rename(curve.id.new = curve.id)
	
	
	


unique(tdata_var$curve.id)
unique(predictions2$curve.id.new)

setdiff(unique(tdata_var$curve.id), unique(predictions2$curve.id.new))

write_csv(all_var, "data-processed/all_var.csv")
write_csv(tdata_var, "data-processed/tdata_var.csv")
write_csv(predictions2, "data-processed/predictions2.csv")
write_csv(data_sel2, "data-processed/data_sel2.csv")

predictions2 %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = mean_temp, y = rate), data = data_sel2, color = "cadetblue") +
	geom_point(aes(x = temperature, y = predicted_rate_var), data = all_var, color = "orange") +
	geom_point(aes(x = temp, y = rate), data = tdata_var, color = "lightblue") +
	facet_wrap( ~ curve.id.new, scales = "free") + ylab("Response") + xlab("Temperature (°C)")
ggsave("figures/all-rate-summation.png", width = 12, height = 10)


unique(predictions2$curve.id.new)

data_sel3 <- data_sel2 %>% 
	rename(temperature = mean_temp) %>% 
	mutate(type = "constant observed") %>% 
	select(curve.id.new, temperature, rate, type)
all_var2 <- all_var %>% 
	rename(rate = predicted_rate_var) %>% 
	mutate(type = "variable predicted") %>% 
	select(curve.id.new, temperature, rate, type)
tdata_var2 <- tdata_var %>% 
	rename(temperature = temp) %>% 
	mutate(type = "variable observed") %>% 
	select(curve.id.new, temperature, rate, type)


all_responses <- bind_rows(data_sel3, all_var2, tdata_var2) 


### trying to fix the mismatch in temperatures for the 168 variable study.
all_responses2 <- all_responses %>% 
	filter(type == "variable observed", curve.id.new == "168_fecundity/body weight_n. nevadensis_116.4236111") %>% View
	mutate(temperature = ifelse(type == "variable observed" &
								curve.id.new == "168_fecundity/body weight_n. nevadensis_116.4236111" &
								temperature == "20", "19.86389", temperature)) %>% 
	mutate(temperature = ifelse(type == "variable observed" &
									curve.id.new == "168_fecundity/body weight_n. nevadensis_116.4236111" &
									temperature == "33", "32.79884", temperature)) 

str(all_responses2)
str(all_responses)

library(beyonce)
colors <- c("cadetblue", "darkgoldenrod1", "darkorange2")

predictions2 %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = temperature, y = rate, color = type), data = all_responses) +
	facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	scale_color_manual(values = colors)
ggsave("figures/all-rate-summation-color.png", width = 14, height = 10)




### get the study 111 data in order
s111 <- filter(all_responses, grepl("Dros",curve.id.new)) %>% 
	filter(type != "variable observed")

thing2 <- filter(tdata2, study_ID == "111", trait == "fecundity", temp_regime == 1) %>% 
	rename(curve.id.new = curve.id) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "_0_Begin Miami flies", "")) %>% 
	mutate(curve.id.new = str_replace(curve.id.new, "_0_Begin New Jersey flies", "")) %>% 
	filter(mean_temp %in% c(18, 25)) %>% 
	mutate(type = "variable observed") %>% 
	mutate(stoch2 = treatment_name_study) %>% 
	mutate(stoch2 = str_replace(treatment_name_study, "stochastic_", "")) %>% 
	filter(stoch2 == mean_temp) %>% 
	rename(temperature = mean_temp) %>% 
	select(curve.id.new, temperature, rate, type)

s111b <- bind_rows(s111, thing2)

predictions2 %>%
	filter(grepl("Dros",curve.id.new)) %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line(color = "cadetblue", size = 1) +
	geom_point(aes(x = temperature, y = rate, color = type),
			   data = s111b, size = 2.5) +
	geom_point(aes(x = temperature, y = rate), shape = 1, color = "black",
		   data = s111b, size = 2.5) +
	# geom_point(aes(x = mean_temp, y = rate), shape = 1, color = "black",
	# 		   data = thing2 )+ 
	facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	# scale_color_brewer(type = "qual", palette = "Dark2") +
	scale_color_manual(values = colors) 
ggsave("figures/study111-predictions.png", width = 12, height = 4)

temps18_111 <- read_csv("Data/study111-temps-18C-mean.csv")
temps25_111 <- read_csv("Data/study111-temps-25C-mean.csv")


all_111_temps <- bind_rows(temps18_111, temps25_111) %>% 
	mutate(study_ID = "111")

all_111_temps %>% 
	ggplot(aes(x = day, y = temperature, color = factor(mean_temperature))) + geom_line() +
	geom_hline(yintercept = 18, color = "grey", linetype = "dashed") +
	geom_hline(yintercept = 25, color = "black", linetype = "dashed") +
	scale_color_manual(values = c("grey", "black"), name = "Mean temperature (°C)") +
	ylab("Temperature (°C)") + xlab("Days")
ggsave("figures/study111-temps.png", width = 12, height = 6)


fluc150 <- read_csv("Data/study150-temperature-fluctuations.csv")

fluc150 %>% 
	ggplot(aes(x = day, y = temperature)) + geom_line() +
	# geom_hline(yintercept = 18, color = "grey", linetype = "dashed") +
	# geom_hline(yintercept = 25, color = "black", linetype = "dashed") +
	scale_color_manual(values = c("grey", "black"), name = "Mean temperature (°C)") +
	ylab("Temperature (°C)") + xlab("Days")

fluc168 <- read_csv("Data/study168-temperature-fluctuations.csv")

fluc168 %>% 
	ggplot(aes(x = hour, y = temperature, color = DTR, group = min_temp)) + geom_line() +
	# geom_hline(yintercept = 18, color = "grey", linetype = "dashed") +
	# geom_hline(yintercept = 25, color = "black", linetype = "dashed") +
	ylab("Temperature (°C)") + xlab("Days")

study168_temps <- read_csv("Data/study168-temperature-fluctuations.csv")


study168_temps %>% 
	rename(fluctuation = flutation_type) %>% 
	mutate(realized_temp = temperature + min_temp) %>% 
	mutate(unique_regime = paste(fluctuation, min_temp, sep = "_")) %>%  
	ggplot(aes(x = hour, y = realized_temp, color = unique_regime, group = unique_regime)) + geom_line()


study168_tempsb <- study168_temps %>% 
	rename(fluctuation = flutation_type) %>% 
	mutate(realized_temp = temperature + min_temp) %>% 
	mutate(unique_regime = paste(fluctuation, min_temp, sep = "_")) %>% 
	mutate(curve.id.list = fits168$curve.id.list[[1]]) 

study168_tempsb %>% 
	ggplot(aes(x = hour, y = realized_temp, color = unique_regime)) + geom_point() +
	geom_line() +
	ylab("Temperature (°C)") + xlab("Hour") +
	scale_color_viridis_d(end = 0.9, name = "Temperature regime")
ggsave("figures/study168-temp-regime.png", width = 12, height = 6)

	

predictions2 %>%
	filter(grepl("nevad",curve.id.new)) %>% 
	distinct()

all_responses %>% 
	filter(grepl("Dros",curve.id.new)) %>% 
	filter(type == "variable observed") %>% View

thing1 <- tdata_var %>% 
	rename(temperature = temp) %>% 
	mutate(type = "variable observed") %>% 
	filter(grepl("Dros",curve.id.new)) %>% 
	filter(type == "variable observed") 



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
	# geom_point(aes(x = temperature, y = growth.rate), data = dat.full, shape = 1, size = 2, color = "grey") +
	# geom_point(aes(x = temperature, y = growth.rate), data = tdata_var, shape = 1, size = 2, color = "cadetblue") +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	# geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate, color = "purple") +
	# geom_line(aes(x = temperature, y = predicted_rate), data = fits_above_zero, color = "cadetblue") +
	xlim(0, 40) + ylab("Rate") + xlab("Temperature (°C)") + ylim(-3, 110) +
	# geom_point(aes(x = mean_temp, y = rate, color = treatment_name_study),
	# 		   data = filter(tdata2, study_ID == "111", trait == "fecundity", temp_regime == 1), size = 2.5)  +
	geom_point(aes(x = temperature, y = rate, color = type),
			   data = filter(all_responses, grepl("Dros",curve.id.new)), size = 2.5) +
	# geom_point(aes(x = temperature, y = growth.rate), data = all_preds, color = "purple") +
	facet_wrap( ~ curve.id.new, scales = "free") 



### now let's bring in the kingsolver data

tdata_king <- read_csv("Data/Kingolver-constant-25.csv") %>% 
	mutate(type = "constant observed") %>% 
	mutate(growth_rate = exp(growth_rate))

king_fits <- read_csv("data-processed/norberg-fits-kingolver.csv") %>% 
	mutate(study_ID = "kingsolver")

king_var <- read_csv("Data/Kingsolver-variable-growth.csv") %>% 
	mutate(temperature = round(mean_temp, digits = 0)) %>% 
	mutate(type = "variable observed") %>% 
	mutate(growth_rate_g = exp(growth_rate_g)) %>% 
	mutate(growth_rate = growth_rate_g)

k_temp_data <- read_csv("Data/Kingolver-variable-temps.csv") %>% 
	rename(mean_temp_25 = temperature) %>% 
	mutate(mean_temp_20 = mean_temp_25 - 5) %>% 
	mutate(mean_temp_30 = mean_temp_25 + 5) %>% 
	gather(key = mean_temp, value = temperature, mean_temp_30, mean_temp_20, mean_temp_25) %>% 
	mutate(study_ID = "kingsolver")




all_king <- left_join(k_temp_data, king_fits) %>% 
	mutate(predicted_growth = a.list*exp(b.list*temperature)*(1-((temperature-z.list)/(w.list/2))^2)) %>% 
	group_by(fluctuation, mean_temp) %>% 
	summarise(growth_rate = mean(predicted_growth)) %>% 
	mutate(type = "variable predicted") %>% 
	mutate(temperature = str_replace(mean_temp, "mean_temp_", "")) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	select(-mean_temp)


fits_split <- king_fits %>% 
	filter(!is.na(topt.list)) %>% 
	split(.$curve.id.list)


fits_king <- fits_split %>% 
	map_df(prediction_function, .id = "curve.id")

fits_above_zero_king <- fits_king %>% 
	filter(predicted_rate > 0, predicted_rate < 150)


king_points <- bind_rows(all_king, king_var, tdata_king)


king_predicted <- read_csv("Data/Kingsolver-variable-growth-predicted-relative.csv") %>% 
	gather(key = type, value = growth_rate, 3:5) %>% 
	filter(type %in% c("observed", "predicted_growth_rate")) %>% 
	mutate(type = str_replace(type, "observed", "variable observed")) %>% 
	mutate(type = str_replace(type, "predicted_growth_rate", "variable predicted"))


king_points_adj <- king_points %>% 
	mutate(growth_rate = ifelse(type == "variable observed", growth_rate + 0.37, growth_rate)) 

fits_above_zero_king %>% 
	# mutate(predicted_rate = predicted_rate - 0.4) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(data = king_points_adj, aes(x = temperature, y = growth_rate, color = type), size =2) +
	geom_point(data = king_points_adj, aes(x = temperature, y = growth_rate), size =2, shape = 1) +
	# geom_point(data = king_predicted, aes(x = temperature, y = growth_rate, color = type)) +
	scale_color_manual(values = colors) +
	ylab("Growth rate (mg/day)") + xlab("Temperature (°C)")
ggsave("figures/kingsolver-tpc-plot.png", width = 6, height = 3.5)
