

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
	filter(temp_regime == 0) %>% 
	ggplot(aes(x = mean_temp_calculated, y = response, color = trait)) + geom_point() +
	facet_wrap(trait ~ study_ID, scales = "free")
ggsave("figures/constant-temp-plots.png", width = 16, height = 8)
