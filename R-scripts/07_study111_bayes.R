
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



