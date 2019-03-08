

## Marta's plotting of the TPCs based on Joey's Code

library(tidyverse)
library(cowplot)

tdata <- read_csv("Data/ExtractedDataAllStudies.csv")

tdata2 <- tdata %>%
	mutate(study_digits = nchar(study_ID)) %>%
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>%
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>%
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) %>%
	filter(study_digits == 2)

# tdata2 %>% 
# 	filter(temp_regime == 0, study_ID %in% c(111)) %>% 
# 	ggplot(aes(x = mean_temp, y = response)) + 
# 	geom_point(size = 2, color = "cadetblue") +
# 	geom_point(size = 2, shape = 1, color = "black") +
# 	facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
# ggsave("Figures/constant-temp-plots-subset-m.png", width = 4, height = 3)
# 
# tdata2 %>% 
# 	filter(study_ID == 111) %>% 
# 	ggplot(aes(x = mean_temp_calculated, y = response, color = factor(temp_regime))) + geom_point() +
# 	facet_wrap(study_ID ~ trait, ncol = 3, scales = "free") + xlab("Temperature (°C)")

unique(tdata2$study_ID)
tdata.58 <- subset(tdata2, study_ID == 58)
tdata.61 <- subset(tdata2, study_ID == 61)
tdata.67 <- subset(tdata2, study_ID == 67)
tdata.90 <- subset(tdata2, study_ID == 90)
tdata.91 <- subset(tdata2, study_ID == 91)

unique(tdata.58$trait)
unique(tdata.61$trait)
unique(tdata.67$trait)
unique(tdata.90$trait)
unique(tdata.91$trait)

tdata.58 %>%
  filter(temp_regime == 0) %>% 
  ggplot(aes(x = mean_temp, y = response)) +
	geom_point(size = 2, color = "cadetblue") +
	geom_point(size = 2, shape = 1, color = "black") +
	facet_wrap(study_ID ~ trait, scales = "free", ncol = 4) + xlab("Temperature (°C)") +ylab("Response")
ggsave("Figures/constant-temp-plots-58.png", width = 20, height = 25)

tdata.61 %>%
  filter(temp_regime == 0) %>% 
  ggplot(aes(x = mean_temp, y = response)) +
  geom_point(size = 2, color = "cadetblue") +
  geom_point(size = 2, shape = 1, color = "black") +
  facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
ggsave("Figures/constant-temp-plots-61.png", width = 8, height = 6)

tdata.67 %>%
  filter(temp_regime == 0) %>% 
  ggplot(aes(x = mean_temp, y = response)) +
  geom_point(size = 2, color = "cadetblue") +
  geom_point(size = 2, shape = 1, color = "black") +
  facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
ggsave("Figures/constant-temp-plots-67.png", width = 12, height = 12)

tdata.90 %>%
  filter(temp_regime == 0) %>% 
  ggplot(aes(x = mean_temp, y = response)) +
  geom_point(size = 2, color = "cadetblue") +
  geom_point(size = 2, shape = 1, color = "black") +
  facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
ggsave("Figures/constant-temp-plots-90.png", width = 4, height = 4)

tdata.91 %>%
  filter(temp_regime == 0) %>% 
  ggplot(aes(x = mean_temp, y = response)) +
  geom_point(size = 2, color = "cadetblue") +
  geom_point(size = 2, shape = 1, color = "black") +
  facet_wrap(study_ID ~ trait, scales = "free", ncol = 3) + xlab("Temperature (°C)") +ylab("Response")
ggsave("Figures/constant-temp-plots-91.png", width = 10, height = 8)
