

library(tidyverse)

colors <- c("cadetblue", "darkgoldenrod1", "darkorange2")

unique(data_sel2$curve.id)

data_sel <- read_csv("data-processed/tdata-sel.csv")
pred_c <- read_csv("data-processed/pred_c.csv") 

all_var <- read_csv("data-processed/all_var.csv")
tdata_var <- read_csv("data-processed/tdata_var.csv")


data_sel2 <- data_sel %>% 
	filter(curve.id %in% c(pred_c$curve.id.new)) %>% 
	mutate(curve.id.new = curve.id)

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


predictions2 <- read_csv("data-processed/predictions2.csv")


## study 168
predictions2 %>%
	filter(grepl("168", curve.id.new)) %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = temperature, y = rate, color = type), data = filter(all_responses, grepl("168", curve.id.new))) +
	geom_point(aes(x = temperature, y = rate), data = filter(all_responses, grepl("168", curve.id.new)), shape = 1) +
	# facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	scale_color_manual(values = colors) + ylab("Fecundity/body weight")
ggsave("figures/study168.png", width = 6, height = 3)

## study 119
predictions2 %>%
	filter(grepl("119", curve.id.new)) %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = temperature, y = rate, color = type), data = filter(all_responses, grepl("119", curve.id.new))) +
	geom_point(aes(x = temperature, y = rate), data = filter(all_responses, grepl("119", curve.id.new)), shape = 1) +
	facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	scale_color_manual(values = colors) + ylab("Weight gain")
ggsave("figures/study119.png", width = 8, height = 3)


## study 150
predictions2 %>%
	filter(grepl("150", curve.id.new)) %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = temperature, y = rate, color = type), data = filter(all_responses, grepl("150", curve.id.new))) +
	geom_point(aes(x = temperature, y = rate), data = filter(all_responses, grepl("150", curve.id.new)), shape = 1) +
	facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	scale_color_manual(values = colors) + ylab("Duration of incubation")
ggsave("figures/study150.png", width = 6, height = 3)


## study 159
unique(predictions2$curve.id.new)

	predictions2 %>%
	filter(grepl("159", curve.id.new)) %>% 
	filter(predicted_rate > 0, predicted_rate < 150) %>% 
	ggplot(aes(x = temperature, y = predicted_rate)) + geom_line() +
	geom_point(aes(x = temperature, y = rate, color = type), data = filter(all_responses, grepl("159", curve.id.new))) +
	geom_point(aes(x = temperature, y = rate), data = filter(all_responses, grepl("159", curve.id.new)), shape = 1) +
	facet_wrap( ~ curve.id.new, scales = "free") +
	ylab("Response") + xlab("Temperature (°C)") +
	scale_color_manual(values = colors) + ylab("Body size")
ggsave("figures/study159.png", width = 12, height = 7)




