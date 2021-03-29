

### Revisiting the analysis in March 2021 (eeks)
### first need to merge the latest extracted data file with Maggie's data and Marta's data and the kingsolver data
### then go through and make sure have merged all the extracted data
### the go back and see how many papers still need to be extracted
### then come up with a gameplan for how to analyze them

library(readxl)
library(tidyverse)
library(janitor)


edata <- read.csv("Data/ExtractedDataAllStudies.csv") %>% 
	clean_names() %>% 
	mutate(temp_range = as.character(temp_range)) %>% 
	mutate(dtr = as.character(dtr)) %>% 
	mutate(min_time = as.character(min_time)) %>% 
	mutate(max_time = as.character(max_time)) %>% 
	mutate(sample_size = as.character(sample_size)) 

unique(edata$variance_type)
mdata <- read_excel("Data/vector-bite-extracted-maggie.xlsx") %>% 
	clean_names() %>% 
	mutate(study_id = as.character(study_id)) %>% 
	mutate(min_time = as.character(min_time)) %>% 
	mutate(max_time = as.character(max_time)) %>% 
	mutate(sample_size = as.character(sample_size)) 

View(edata)
names(mdata)
names(edata)

adata <- bind_rows(edata, mdata)
