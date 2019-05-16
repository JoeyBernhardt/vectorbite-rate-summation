

library(tidyverse)
library(temperatureresponse)
library(cowplot)

tdata <- read_csv("Data/ExtractedDataAllStudies.csv")


tdata2 <- tdata %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(study_digits == 3)




tdata2 %>% 
	filter(study_ID == 159, temp_regime == 0) %>% 
	ggplot(aes(x = mean_temp_calculated, y = response, color = factor(temp_regime))) + geom_point() +
	facet_wrap(study_ID ~ trait, ncol = 6, scales = "free") + xlab("Temperature (°C)") + geom_smooth()

snippet <- tdata2 %>% 
	filter(study_ID == 159, temp_regime == 0, trait == "winglength") %>% 
	rename(temp = mean_temp_calculated,
		   rate = response)
snippet_var <- tdata2 %>% 
	filter(study_ID == 159, temp_regime == 1, trait == "winglength") %>% 
	rename(temp = mean_temp_calculated,
		   rate = response) %>% 
	select(temp, rate, everything())

equ8

output <- with(snippet, fitmodellist(temp=temp, rate=rate))
output_var <- with(snippet_var, fitmodellist(temp=temp, rate=rate))

out <- with(snippet, equ10(temp=temp, rate=rate))

eq8_fit <- output %>% 
	filter(model == "equ08")

eq8_fit_var <- output_var %>% 
	filter(model == "equ08")

eq4 <-  function(t, R = 0.001987) {
	t <- t + 273.15
	y <- eq4_fit$estimate[eq4_fit$term == "a"] * exp(-eq4_fit$estimate[eq4_fit$term == "b"]/(R * t))
	- eq4_fit$estimate[eq4_fit$term == "c"] * exp(-eq4_fit$estimate[eq4_fit$term == "d"]/(R * t))
}

eq8 <-  function(t) {
	tref <- eq8_fit$estimate[eq8_fit$term == "tref"] 
	y <- eq8_fit$estimate[eq8_fit$term == "a"] * exp(-0.5 * ((t - tref)/eq8_fit$estimate[eq8_fit$term == "b"])^2)
}

eq8_var <-  function(t) {
	tref <- eq8_fit_var$estimate[eq8_fit_var$term == "tref"] 
	y <- eq8_fit_var$estimate[eq8_fit_var$term == "a"] * exp(-0.5 * ((t - tref)/eq8_fit_var$estimate[eq8_fit$term == "b"])^2)
}

exp(-10)

eq4_fit$estimate[eq4_fit$term == "CTmin"]

equ8

temps <- seq(0, 40, length = 50)

predictions <- sapply(temps, eq8)
predicted_rate <- data.frame(temperature = temps, predicted_rate = predictions)


predictions_var <- sapply(temps, eq8_var)
predicted_rate_var <- data.frame(temperature = temps, predicted_rate = predictions_var)


ggplot(aes(x = temp, y = rate), data = snippet) + geom_point(shape = 1, size = 2, color = "grey") +
	geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate) +
	geom_line(aes(x = temperature, y = predicted_rate), data = predicted_rate_var, color = "purple") +
	geom_vline(xintercept =  eq8_fit$estimate[eq8_fit$term == "topt"], linetype = 2, color = "grey")  +
	geom_point(aes(x = temp, y = rate), data = snippet_var, color = "purple") + ylab("Rate") + xlab("Temperature (°C)")
