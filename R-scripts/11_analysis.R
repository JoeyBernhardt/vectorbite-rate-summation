

### Revisiting the analysis in March 2021 (eeks)
### first need to merge the latest extracted data file with Maggie's data and Marta's data and the kingsolver data
### then go through and make sure have merged all the extracted data
### the go back and see how many papers still need to be extracted
### then come up with a gameplan for how to analyze them

library(readxl)
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())


edata <- read.csv("Data/ExtractedDataAllStudies.csv") %>% 
	clean_names() %>% 
	mutate(temp_range = as.character(temp_range)) %>% 
	mutate(dtr = as.character(dtr)) %>% 
	mutate(min_time = as.character(min_time)) %>% 
	mutate(max_time = as.character(max_time)) %>% 
	mutate(sample_size = as.character(sample_size)) %>% 
	select(1:9, trait, units, response) %>% 
	mutate(dataset = "vb_team_data")

unique(edata$variance_type)
mdata <- read_excel("Data/vector-bite-extracted-maggie.xlsx") %>% 
	clean_names() %>% 
	mutate(study_id = as.character(study_id)) %>% 
	mutate(min_time = as.character(min_time)) %>% 
	mutate(max_time = as.character(max_time)) %>% 
	mutate(sample_size = as.character(sample_size)) %>% 
	select(1:9, trait, units, response) %>% 
	mutate(dataset = "maggie_data")

tdata_king <- read_csv("Data/Kingolver-constant-25.csv") %>% 
	mutate(type = "constant observed") %>% 
	mutate(growth_rate = exp(growth_rate)) %>% 
	mutate(temp_regime = 0)

king_var <- read_csv("Data/Kingsolver-variable-growth.csv") %>% 
	mutate(temperature = round(mean_temp, digits = 0)) %>% 
	mutate(type = "variable observed") %>% 
	mutate(growth_rate_g = exp(growth_rate_g)) %>% 
	mutate(growth_rate = growth_rate_g) %>% 
	mutate(temp_regime = 1)

all_king <- bind_rows(tdata_king, king_var) %>% 
	select(-mean_temp) %>% 
	rename(mean_temp = temperature) %>% 
	rename(response = growth_rate) %>% 
	mutate(study_id = "kingsolver") %>% 
	mutate(trait = "growth rate")

View(edata)
names(mdata)
names(edata)

adata <- bind_rows(edata, mdata, all_king) %>% 
	mutate(mean_temp_calculated = (min_temp + max_temp)/2) %>% 
	mutate(mean_temp_calculated = ifelse(temp_regime == 0, mean_temp, mean_temp_calculated)) %>% 
	mutate(mean_temp_calculated = ifelse(is.na(mean_temp_calculated), mean_temp, mean_temp_calculated)) 

View(adata)

unique(adata$study_id)


adata2 <- adata %>% 
  
  filter(!is.na(response)) %>% 
  
  dplyr::group_by(study_id, trait) %>% 
  mutate(max_response = max(response)) %>% 
  mutate(scaled_response = response/max_response) %>% 
  mutate(unique_experiment = paste(study_id, trait, sep = "_"))


plot1 <- adata2 %>% 
  ggplot(aes(x = mean_temp_calculated, y = scaled_response, group = unique_experiment, color = factor(temp_regime))) +
  geom_point(size = 2) +
  facet_wrap(study_id ~ trait, ncol = 6, scales = "free")+
  xlab("Temperature (°C)")
plot1
ggsave(plot = plot1, "figures/all_studies_plots_LC.png", width = 35, height = 45)


num_points <- adata2 %>% 
	dplyr::group_by(study_id, trait, temp_regime) %>% 
	tally() %>% 
	mutate(keep = ifelse(n > 4, "keep", "drop"))

number_of_points <- adata2 %>% 
	dplyr::group_by(study_id, trait, temp_regime) %>% 
	tally() %>%
	mutate(treatment = ifelse(temp_regime == 0, "constant", "fluctuating"))

write_csv(number_of_points, "data-processed/number_of_points.csv")
	

adata3 <- left_join(adata2, num_points) %>% 
	filter(keep == "keep")

unique(adata2$trait)

#write_csv(adata3, "data-processed/alldata_scaled.csv")
