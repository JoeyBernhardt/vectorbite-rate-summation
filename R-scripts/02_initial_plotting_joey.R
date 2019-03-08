

## Joey's plotting of the TPCs


library(tidyverse)
library(cowplot)


tdata <- read_csv("Data/ExtractedDataAllStudies.csv")


tdata2 <- tdata %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>% 
	filter(study_digits == 3)


tdata2 %>% 
	filter(temp_regime == 0, study_ID %in% c(111)) %>% 
	ggplot(aes(x = mean_temp_calculated, y = response)) + 
	geom_point(size = 2, color = "cadetblue") +
	geom_point(size = 2, shape = 1, color = "black") +
	facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
ggsave("figures/constant-temp-plots-subset-111.png", width = 4, height = 3)

tdata2 %>% 
	filter(study_ID == 111) %>% 
	ggplot(aes(x = mean_temp_calculated, y = response, color = factor(temp_regime))) + geom_point() +
	facet_wrap(study_ID ~ trait, ncol = 3, scales = "free") + xlab("Temperature (°C)")
