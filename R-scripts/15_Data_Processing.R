## Lilian Chan, University of Guelph
## Rate summation projection
##
## Purpose: merge data from the vb_team, maggie's data, and data from Slein et 
## al. 2023
## 
##
## Table of content:
##    0. Set-up workspace
##
##    1. Organizing and combining datasets
##        A. vb_team
##        B. Maggie's data (including Kingsolver paper)
##        C. Slein et al. 2023
##        D. Combine datasets
##
##    2. Filter suitable papers/data
##        A. Plot the raw data for each curve id
##        B. Filter papers with >= 3 constant treatment
##        C. Data processing


##########
###### 0. Set-up workspace ----
##########
library(tidyverse)
library(readxl)
library(janitor)
library(cowplot)

setwd("~/Documents/UofG/vectorbite-rate-summation")


##########
###### 1A. Organizing vb_team data ----
##########

edata <- read.csv("Data/ExtractedDataAllStudies.csv") %>% 
  clean_names() %>% 
  mutate(temp_range = as.character(temp_range)) %>% 
  mutate(dtr = as.character(dtr)) %>% 
  mutate(min_time = as.character(min_time)) %>% 
  mutate(max_time = as.character(max_time)) %>% 
  mutate(dtr = as.numeric(dtr)) %>% 
  mutate(sample_size = as.character(sample_size)) %>% 
  mutate(variance_value = as.numeric(variance_value)) %>% 
  mutate(qual_measure = as.factor(qual_measure)) %>% 
  mutate(ind_pop = as.factor(ind_pop)) %>% 
  mutate(ls_treat = as.factor(ls_treat)) %>% 
  mutate(lab_wild = as.factor(lab_wild)) %>% 
  mutate(dataset = "vb_team_data")

colnames(edata)

##########
###### 1B. Organizing Maggie's data ----
##########

mdata <- read_excel("Data/vector-bite-extracted-maggie.xlsx") %>% 
  clean_names() %>% 
  mutate(study_id = as.character(study_id)) %>% 
  mutate(temp_init = as.character(temp_init)) %>% 
  mutate(min_time = as.character(min_time)) %>% 
  mutate(max_time = as.character(max_time)) %>% 
  mutate(dtr = as.numeric(dtr)) %>% 
  mutate(sample_size = as.character(sample_size)) %>% 
  mutate(variance_type = as.character(variance_type)) %>% 
  mutate(qual_measure = as.factor(qual_measure)) %>% 
  mutate(ind_pop = as.factor(ind_pop)) %>% 
  mutate(ls_treat = as.factor(ls_treat)) %>% 
  mutate(lab_wild = as.factor(lab_wild)) %>% 
  mutate(dataset = "maggie_data")

colnames(mdata)

## Kingsolver paper (this is Kingsolver et al. 2015/ study id 16)
king_var <- read_csv("Data/Kingsolver-variable-growth.csv") %>% 
  clean_names() %>% 
  mutate(study_id = "kingsolver") %>% 
  mutate(temp_regime = 1) %>% 
  mutate(mean_temp = round(mean_temp, digits = 0)) %>% 
  mutate(min_temp = mean_temp - fluctuation) %>% 
  mutate(max_temp = mean_temp + fluctuation) %>% 
  rename(dtr = fluctuation) %>% 
  mutate(fluc_type = 2) %>% 
  mutate(ramp_type = 1) %>% 
  mutate(trait = "growth rate") %>% 
  mutate(units = "g/day") %>% 
  mutate(response = exp(growth_rate_g)) %>% 
  mutate(qual_measure = as.factor(1)) %>% 
  mutate(genus = "Manduca") %>% 
  mutate(species = "sexta") %>% 
  mutate(dataset = "maggie_data")

## Reorder the columns
king_var <- king_var %>% 
  select("study_id", "temp_regime", "mean_temp", "min_temp", "max_temp", "dtr", 
         "fluc_type", "ramp_type", "trait", "units", "response", 
         "qual_measure", "genus", "species", "dataset")


##########
###### 1C. Organizing Slein et al. 2023's data ----
##########

## get the study_id of the usable papers from Slein et al. 2023
# slein_needed <- read_excel("Data Papers/study_info.xlsx", sheet = "Sheet1") %>% 
#   filter(dataset == "Slein_cited_ref" | dataset == "Slein_original")
# 
# study_id_needed <- c(slein_needed$study_id)

## read the data from Slein et al. 2023
slein_data <- read.csv("Data/mslein/data/rawdata23dec22.csv") %>% 
  clean_names() #%>% 
  #filter(is.element(study_id, study_id_needed))

unique(slein_data$study_id)


## select columns that are equivalent to the columns in the original dataset
sdata <- slein_data %>% 
  select(study_id, # study_id
         secondary_temp, # temp_init
         mean_temp_constant, # mean_temp
         flux_range, # dtr
         flux_pattern, # fluc_type
         period_flux, #mean_time, min_time, max_time
         resp_def, # trait
         resp_units, # units
         constant_resp, # response
         flux_resp, # response
         constant_samp, # sample_size (need double check)
         flux_samp, # sample_size (need double check)
         variance_type, # variance_type
         constant_variance, # variance_value
         flux_variance, # variance_value
         resp_quality, # qual_measure
         org_level, # ind_pop (slein has community(2) and ecosystem (3) levels)
         genus, # genus
         species,# species
         exp_age, # ls_treat (need double check)
         size # body_size (DIFFERENT CODES!!!!! DO NOT MERGE)
         ) 


## reshape the dataset such that each row is a measurement
sdata <- sdata %>% 
  pivot_longer(cols = c(constant_resp, flux_resp, 
                        constant_samp, flux_samp, 
                        constant_variance, flux_variance
                        ),
               names_to = c("dtr1", ".value"),
               names_sep = "_")


## Change DTR into 0 (if constant) or flux_range (if fluctuate)
sdata <- sdata %>% 
  mutate(dtr = ifelse(dtr1 == "constant", 0, flux_range))


## Rename the columns
colnames(sdata) <- c("study_id", "temp_init", "mean_temp", "flux_range", 
                     "fluc_type", "period_flux", "trait", "units", "variance_type",
                     "qual_measure", "ind_pop", "genus", "species", "ls_treat", 
                     "body_size", "dtr1", "response", "sample_size", 
                     "variance_value", "dtr")

## Reorder the columns
sdata <- sdata %>% 
  select("study_id", "temp_init", "mean_temp", "dtr", 
         "fluc_type",# "period_flux",
         "trait", "units", "response", 
         "sample_size", "variance_type", "variance_value", "qual_measure", 
         "ind_pop", "genus", "species", "ls_treat", "body_size") %>% 
  mutate(temp_init = as.character(temp_init)) %>% 
  mutate(sample_size = as.character(sample_size)) %>% 
  mutate(variance_type = as.character(variance_type)) %>% 
  mutate(qual_measure = as.factor(qual_measure)) %>% 
  mutate(ind_pop = as.factor(ind_pop)) %>% 
  mutate(ls_treat = as.factor(ls_treat)) %>% 
  mutate(dataset = "slein_data")


##########
###### 1D. Combine datasets ----
##########

## Combine all data into one dataset
adata <- bind_rows(edata, mdata, king_var, sdata) %>% 
  relocate(additional_complexity, .after = species) %>% 
  # if constant treatment -> DTR = 0
  mutate(dtr_cal = ifelse(temp_regime == 0, 0, max_temp - min_temp)) %>% 
  # if we couldn't calculate DTR from min and max temp but DTR is provided, use the original DTR
  mutate(dtr_cal = ifelse((is.na(dtr_cal) & !is.na(dtr)), dtr, dtr_cal)) %>% 
  relocate(dtr_cal, .after = dtr)

View(adata)

colnames(adata)

problematic <- adata %>% 
  filter(dtr != dtr_cal)

unique(problematic$study_id) # Need to double check these papers

unique(adata$study_id) # There are 109 papers in this dataset

## Add a curve_id that is unique for each study_id, trait, species, and additional complexity
adata <- adata %>% 
  group_by(study_id, trait, species, additional_complexity) %>% 
  mutate(curve_id = cur_group_id()) %>% 
  ungroup() %>% 
  relocate(curve_id, .before = study_id) %>% 
  arrange(curve_id)

## Save the dataset
#write_csv(adata, "data-processed/all_data_20250430.csv")


##########
###### 2A. Plot the raw data for each curve id ----
##########

data <- read_csv("data-processed/all_data_20250430.csv")

# data$dtr <- as.factor(data$dtr)
# data$dtr_cal <- as.factor(data$dtr_cal)


# Use dtr_cal as DTR
data <- data[, -9] %>% 
  rename(dtr = dtr_cal)

# # Create a directory for output files
# output_dir <- "figures/data-paper-plots"
# 
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# 
# 
# # Generate multiple plots and save them
# unique_ids <- unique(data$curve_id)  # Store unique IDs for iteration
# 
# for (i in seq_along(unique_ids)) {
#   # Create subset of data
#   mini_df <- data %>% filter(curve_id == unique_ids[i])
#   
#   # Create a ggplot
#   p <- ggplot(mini_df, aes(x = mean_temp, y = response, colour = dtr)) +
#     geom_point() +
#     labs(title = paste("Curve ID", unique(mini_df$curve_id), "Study ID",
#                        unique(mini_df$study_id), unique(mini_df$trait))) +
#     theme_bw()
#   
#   # Define file path
#   file_name <- paste0(output_dir, "/plot_", unique(mini_df$curve_id), ".png")
#   
#   # Save the plot
#   ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300)
#   
#   print(paste("Saved:", file_name))
# }

##########
###### 2B. Filter papers with >= 3 constant treatment ----
##########

## Since the Slein dataset is incomplete, I will focus on maggie_data and vb_team_data for now

data <- data %>%
  filter(dataset != "slein_data")

# num_of_const <- data %>% 
#   filter(dtr == 0) %>% 
#   group_by(curve_id) %>% 
#   tally()
# 
# ## Remove papers with <3 data from constant temp treatment ----
# curve_needed <- num_of_const %>% 
#   filter(n >= 3)
# 
# curve_needed <- c(curve_needed$curve_id) # 432 curves are suitable (108 from maggie and vb_team)
# 
# curve_needed
# 
# ## Filter dataset
# data_filter <- data %>% 
#   filter(is.element(curve_id, curve_needed)) %>% 
#   arrange(curve_id) %>% 
#   select(curve_id, study_id, temp_regime, mean_temp, dtr, response, trait, units,
#          additional_complexity, ls_treat, ls_measure)

unique(data$study_id) # from 79 papers (25 from maggie and vb_team)

# view(data_filter %>% filter(ls_treat != ls_measure))

unique(data$trait)


## Exclude study ID 11, 61 (no constant treatment); 111, 58 (acclimation not acute); 
## 12 (no DTR); 119 (only keep water lettuce data curve id 13); 124 (see notes in study_info); 
## 16 (problems with digitizing data), kingsolver (repeat of 16)
data <- data[-(which(data$study_id == "11")), ]
data <- data[-(which(data$study_id == "61")), ]
data <- data[-(which(data$study_id == "111")), ]
data <- data[-(which(data$study_id == "58")), ]
data <- data[-(which(data$study_id == "12")), ]
data <- data[-(which(data$study_id == "119" & data$curve_id != 13)), ]
data <- data[-(which(data$study_id == "124")), ]
data <- data[-(which(data$study_id == "16")), ]
data <- data[-(which(data$study_id == "kingsolver")), ]

## Process data ----

## Study ID 10 (curve_id = 1 & 2) ----
## i. convert to development rate by (1/development time)
### convert those that did not develop into NA
data$response[data$study_id == "10" & data$response == 0] <- NA
data$response[data$study_id == "10"] <- 1/data$response[data$study_id == "10"]

### change trait to development rate and unit to 1/days
data$trait[data$study_id == "10"] <- "development rate"
data$trait_def[data$study_id == "10"] <- "1/time to hatch"
data$units[data$study_id == "10"] <- "1/days"

### Check the data
studyid10 <- data[data$study_id == "10",]

studyid10 %>% ggplot(aes(x = mean_temp, y = response, colour = as.factor(dtr))) +
  geom_point() +
  facet_wrap(~additional_complexity) +
  theme_bw()


## Study ID 13 (curve_id 19-23)----
## According to Fig. 1 caption, mean temperatures of the fluctuating regimes were 22 and 26°C for 18-28°C and 18-34°C treatments respectively.
## Edit the mean_temp for 18-28°C
data$mean_temp[data$study_id == "13" & data$dtr == "10"] <- 22


## Study ID 15 (curve_id 24-25) ----
## This is actually Siddiqui and Barlow et al. 1972
## i. add trait name to 26ºC and 27.5ºC (same as others)
data$trait[data$study_id == "15" & is.na(data$trait)] <- "innate capacity for increase (rm)"

## ii. calculate mean_temp
data$mean_temp[data$study_id == "15" & is.na(data$mean_temp)] <- # select data that requires calculation
  rowMeans(data[data$study_id == "15" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T) 

## iii. one of the max_time is not right
data$max_time[which(data$study_id == "15" & data$max_time == "43200")] <- 720


## Study ID 150 (curve_id 26-30) ----
## min_temp = 22.3°C, max_temp = 32.8°C
data$min_temp[data$study_id == "150" & is.na(data$min_temp)] <- 22.3
data$max_temp[data$study_id == "150" & is.na(data$max_temp)] <- 32.8

## DTR
data$dtr[data$study_id == "150" & is.na(data$dtr)] <- 
  data$max_temp[data$study_id == "150" & is.na(data$dtr)] - 
  data$min_temp[data$study_id == "150" & is.na(data$dtr)]

## Convert incubation duration to rate
data$response[data$study_id == "150" & data$trait == "Duration of incubation"] <- 1/data$response[data$study_id == "150" & data$trait == "Duration of incubation"]
data$units[data$study_id == "150" & data$trait == "Duration of incubation"] <- "1/days"
data$trait[data$study_id == "150" & data$trait == "Duration of incubation"] <- "1/Duration of incubation"
data$trait_def[data$study_id == "150" & data$trait == "1/Duration of incubation"] <- "1/number of days from laying to hatching"

## Study ID 159 (curve_id 31-34) ----
## i. separate data into F and M and 2 different populations
data$additional_complexity[data$study_id == "159"] <- data$notes[data$study_id == "159"]

studyid159 <- data[data$study_id == "159",]

data$additional_complexity[data$study_id == "159"] <- c(rep(c(
  # D. melanogaster
  rep("Cordoba - F", 21), rep("Bordeaux - F", 21), rep("Cordoba - M", 21), rep("Bordeaux - M", 21),
  # D. simulans
  rep("Cordoba - F", 19), rep("Bordeaux - F", 19), rep("Cordoba - M", 19), rep("Bordeaux - M", 19)
  ), 
  # 2 traits
  2))


## ii. calculate mean_temp for fluc (just average min_temp & max_temp; see table 3)
data$mean_temp[data$study_id == "159" & is.na(data$mean_temp)] <-
  rowMeans(data[data$study_id == "159" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T) 

## iii. change units of one of the cells into mm
data$units[data$study_id == "159" & is.na(data$units)] <- "mm"

### check data
studyid159 <- data[data$study_id == "159",]

ggplot() +
  geom_point(data = studyid159 %>% filter(trait == "thorax length"),
             aes(x = mean_temp, y = response, color = as.factor(dtr), shape = additional_complexity)) +
  geom_line(data = studyid159 %>% filter(trait == "thorax length" & dtr == 0), aes(x = mean_temp, y = response, group = additional_complexity)) +
  ylim(0.8,1.16) +
  facet_grid(rows = vars(species)) +
  scale_shape_manual(values=c(16, 15, 21,18))+
  theme_bw()


## Study ID 165 (curve_id 38-39) ----
## i. separate data into 2 insect colony (MA and NC) and 2 host plant (potato, horse-nettle)
### they are found in column "notes" and "x
studyid165 <- data[data$study_id == "165",]
studyid165 <- unite(data[data$study_id == "165",], notes, x, 
                    col = "additional_complexity", sep = "-", remove = F) %>% # combine cells into one vales
  relocate(additional_complexity, .after = species)

data[data$study_id == "165",] <- studyid165

### check data
studyid165[studyid165$trait == "development rate",] %>% 
  ggplot(aes(x = mean_temp, y = response, color = as.factor(dtr))) +
  geom_point() +
  geom_line() +
  ylim(0,0.1) +
  facet_wrap(~additional_complexity)
  theme_bw()

## ii. convert min_time and max_time to minutes
data$min_time[which(data$study_id == "165" & data$min_time == "16")] <- 16*60
data$max_time[which(data$study_id == "165" & data$max_time == "8")] <- 8*60


## Study ID 168 (curve_id 40-45) (SKIP THIS FOR NOW) ----
## Need to calculate mean_temp for each DTR (see Fig. 1 and study168-temperature-fluctuations.csv in Data)
studyid168_temp <- read_csv("Data/study168-temperature-fluctuations.csv")

studyid168_temp %>% 
  group_by(min_temp, DTR) %>% 
  summarise(mean_temp_change = mean(temperature)) %>%
  mutate(mean_temp = mean_temp_change + min_temp)

## FOR NOW: (mean_temp is now (min+max)/2 !!!
data$mean_temp[data$study_id == "168" & is.na(data$mean_temp)] <-
  rowMeans(data[data$study_id == "168" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T)


## Study ID 67 (curve_id 84-92) ----
## i. separate data into 2 different strains
data$additional_complexity[data$study_id == "67"] <- rep(c(rep("Thely", 6), rep("Arrheno", 6)), 9)

## ii. since there is no response for all other traits except tibia length and development rate, remove all other data
data <- data[-(which((data$study_id == "67") & !(data$trait  %in% c("tibia length", "development rate")))), ]

### check data
studyid67 <- data[data$study_id == "67",]

studyid67 %>% 
  ggplot(aes(x = mean_temp, y = response, color = as.factor(dtr), shape = additional_complexity)) +
  geom_point() + 
  geom_line() +
  facet_wrap(~trait, scales = "free") +
  theme_bw()


## Study ID 90 (curve_id 98) ----
## i. separate data into infected and uninfected populations (see notes)
data$additional_complexity[data$study_id == "90"] <- c(rep("infected", 10), rep("uninfected", 10))

### check data
studyid90 <- data[data$study_id == "90",]

studyid90 %>% 
  ggplot(aes(x = mean_temp, y = response, color = as.factor(dtr), shape = additional_complexity)) +
  geom_point() + 
  geom_line() +
  theme_bw()

## ii. edit data - min_time max_time do not look right! (see table S1) (SKIP FOR NOW)


## Study ID VB19-1 (curve_id 104) ----
## i. edit mean_temp (same as constant treatment)
data$mean_temp[data$study_id == "VB19-1" & is.na(data$mean_temp)] <-
  rowMeans(data[data$study_id == "VB19-1" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T)

## ii. Convert development time to development rate
data$response[data$study_id == "VB19-1"] <- 1/data$response[data$study_id == "VB19-1"]

### change trait to development rate and unit to 1/days
data$trait[data$study_id == "VB19-1"] <- "development rate"
data$trait_def[data$study_id == "VB19-1"] <- "1/time to development"
data$units[data$study_id == "VB19-1"] <- "1/days"

### Check the data
studyidVB19_1 <- data[data$study_id == "VB19-1",]

studyidVB19_1 %>% ggplot(aes(x = mean_temp, y = response, colour = as.factor(dtr))) +
  geom_point() +
  theme_bw()


## Study ID VB19-2 (curves_id 105-109) ----
## i. calculate mean_temp (16h*max + 8h*min)/24h
data$mean_temp[data$study_id == "VB19-2" & is.na(data$mean_temp)] <-
  (data$min_temp[data$study_id == "VB19-2" & !is.na(data$min_temp)] * 8 +
  data$max_temp[data$study_id == "VB19-2" & !is.na(data$max_temp)] * 16)/24

## ii. convert min_time and max_time to minutes
data$min_time[which(data$study_id == "VB19-2" & data$min_time == "8h")] <- 8*60
data$max_time[which(data$study_id == "VB19-2" & data$max_time == "16h")] <- 16*60


## Study ID VB19-3 (curve_id 110-111) ----
## i. calculate mean_temp
data$mean_temp[data$study_id == "VB19-3" & is.na(data$mean_temp)] <-
  rowMeans(data[data$study_id == "VB19-3" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T)

## ii. convert development time to development rate
data$response[data$study_id == "VB19-3"] <- 1/data$response[data$study_id == "VB19-3"]

### change trait to development rate and unit to 1/days
data$additional_complexity[data$study_id == "VB19-3"] <- c(rep("female", 21), rep("male", 21))
data$trait[data$study_id == "VB19-3"] <- "development rate"
data$trait_def[data$study_id == "VB19-3"] <- "1/time to development"
data$units[data$study_id == "VB19-3"] <- "1/days"

## iii. convert min_time and max_time to mins
data$min_time[which(data$study_id == "VB19-3" & data$min_time == "12h")] <- 12*60
data$max_time[which(data$study_id == "VB19-3" & data$max_time == "12h")] <- 12*60


## Study ID VB19-4 (curve_id 112) ----
## i. convert min_time and max_time to mins
data$min_time[which(data$study_id == "VB19-4" & data$min_time == "8h")] <- 8*60
data$max_time[which(data$study_id == "VB19-4" & data$max_time == "16h")] <- 16*60

## ii. change the dtr of the fluc treatment
data$dtr[which(data$study_id == "VB19-4")][9] <- 
  data$max_temp[which(data$study_id == "VB19-4")][9] - 
  data$min_temp[which(data$study_id == "VB19-4")][9]

## also the temp_regime
data$temp_regime[which(data$study_id == "VB19-4")][9] <- 1

## iii. convert development time to development rate
data$response[data$study_id == "VB19-4"] <- 1/data$response[data$study_id == "VB19-4"]

### change trait to development rate and unit to 1/days
data$trait[data$study_id == "VB19-4"] <- "development rate"
data$trait_def[data$study_id == "VB19-4"] <- "1/total development time - egg to adult"
data$units[data$study_id == "VB19-4"] <- "1/days"


## Study ID VB19-5 (curve_id 113-117):
## i. calculate mean_temp ((min_temp + max_temp)/2)
data$mean_temp[data$study_id == "VB19-5" & is.na(data$mean_temp)] <-
  rowMeans(data[data$study_id == "VB19-5" & is.na(data$mean_temp), c("min_temp", "max_temp")], na.rm = T)

## ii. convert development time to development rate
data$response[data$study_id == "VB19-5"][23:44] <- 1/data$response[data$study_id == "VB19-5"][23:44]

data$additional_complexity[data$study_id == "VB19-5"] <- c(rep("female", 11), rep("male", 11), rep("female", 11), rep("male", 11), rep(NA, 11))
data$trait[data$study_id == "VB19-5"][23:44] <- "development rate"
data$trait_def[data$study_id == "VB19-5"][23:33] <- "1/time to development"
data$trait_def[data$study_id == "VB19-5"][34:44] <- "1/duration of 50% hatching to enclusion"
data$units[data$study_id == "VB19-5"][23:44] <- "1/days"


## Give a unique curve_id (unique for each study_id, trait, species, and additional complexity)
data <- data %>% 
  group_by(study_id, trait, species, additional_complexity) %>% 
  mutate(curve_id = cur_group_id()) %>% 
  ungroup() %>% 
  relocate(curve_id, .before = study_id) %>% 
  arrange(curve_id) %>% 
  filter(!is.na(response))


## Check the data
num_of_const <- data %>% 
  filter(dtr == 0) %>% 
  group_by(curve_id) %>% 
  summarise(n = n()) %>% 
  arrange(n)

## Remove papers with <4 data from constant temp treatment ----
curve_needed <- num_of_const %>% 
  filter(n >= 4)

curve_needed <- c(curve_needed$curve_id) # 60 from maggie and vb_team
curve_needed

data <- data %>% filter(curve_id %in% curve_needed)
unique(data$curve_id)


## Do the curve_id again (unique for each study_id, trait, species, and additional complexity)
data <- data %>% 
  group_by(study_id, trait, species, additional_complexity) %>% 
  mutate(curve_id = cur_group_id()) %>% 
  ungroup() %>% 
  relocate(curve_id, .before = study_id) %>% 
  arrange(curve_id)

unique(data$curve_id) # 60 curves from vb_team and maggie's data


## Save the dataset
# write_csv(data, "data-processed/vb_maggie_data_20250611.csv")

data <- read_csv("data-processed/vb_maggie_data_20250611.csv")


## Plot each curve again
# Create a directory for output files
output_dir <- "figures/Lilian/data-paper-plots-20250618"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Generate multiple plots and save them
unique_ids <- unique(data$curve_id)  # Store unique IDs for iteration

for (i in seq_along(unique_ids)) {
  # Create subset of data
  mini_df <- data %>% filter(curve_id == unique_ids[i])

  # Create a ggplot
  p <- ggplot(mini_df, aes(x = mean_temp, y = response, colour = as.factor(dtr))) +
    geom_point() +
    labs(title = paste("Curve ID", unique(mini_df$curve_id), "Study ID",
                       unique(mini_df$study_id), unique(mini_df$trait), 
                       unique(mini_df$additional_complexity))) +
    theme_bw()

  # Define file path
  file_name <- paste0(output_dir, "/plot_", unique(mini_df$curve_id), "_", 
                      unique(mini_df$study_id), ".png")

  # Save the plot
  ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300)

  print(paste("Saved:", file_name))
}


## Save a file show the curve_id, study_id, genus-species, additional complexity, trait, trait_def, units
study_info_curves <- data %>% 
  unite(., genus, species, col = "genus_species", sep = " ") %>% 
  select(curve_id, study_id, genus_species, additional_complexity, trait, trait_def, units)

## Each row is each curve_id
study_info_curves <- study_info_curves[!duplicated(study_info_curves),]

# write_csv(study_info_curves, "Data Papers/study_info_curves.csv")
