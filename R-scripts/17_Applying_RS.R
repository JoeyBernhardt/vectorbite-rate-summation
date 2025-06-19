## Lilian Chan, University of Guelph
## Rate summation projection
##
## Purpose: Fit TPCs (Briere, quadratic, Sharpe-Schoolfield (full), Lactin) 
## to the data from constant temperature treatment
## 
##
## Table of content:
##    0. Set-up workspace
##    1. Preparation
##        A. Extract TPC under constant temperature
##        B. Time series for the fluctuating regime
##    2. Rate summation calculations
##    3. Plotting


##########
###### 0. Set-up workspace ----
##########
## load packages
## Load Packages
library(readxl)
library(tidyverse)
library(janitor)
library(progress)
library(ggsci) # colour palettes
library(ggpubr)
library(ggrepel)

setwd("~/Documents/UofG/vectorbite-rate-summation")

## Load functions
source("R-scripts/00_Functions.R")



##########
###### 1A. Extract TPC under constant temperature  ----
##########

## Load data
data <- read_csv("data-processed/vb_maggie_data_20250611.csv")

## Remove Study ID 150 since it is not diurnal fluctuation
data <- data %>% filter(study_id != "150")

## Load the constant temperature TPC of each curve
load("R-scripts/model-output/d_preds_best.Rdata")

# Round the temperature to the nearest 0.1 (to avoid errors when doing joining later)
d_preds_best$temp <- round(d_preds_best$temp, 1)

## temperature gradient to apply rate summation 
temp_grad <- round(seq(0, 50, 0.1), 1) 


## Study IDs:
study_ids <- unique(data$study_id)

## We have data from 15 studies: 10,119,13,15,150,159,165,168,9,90,VB19-1,VB19-2,VB19-3,VB19-4,VB19-5
study_ids

dtr <- data %>% 
	filter(dtr != 0) %>% 
	group_by(study_id) %>% 
	arrange(study_id, dtr) %>%
	distinct(dtr)

## Round dtr to the nearest integer for convenience
dtr$dtr <- round(dtr$dtr,0)

# ## The predictions of the Briere model becomes NaN above Tmax -> change those predictions to 0
# d_preds_best$.fitted[d_preds_best$model_name == "briere" & is.nan(d_preds_best$.fitted)] <- 0


##########
###### 1B. Time series for the thermal fluctuation regime ----
##########

## Generate the temperature profile for a range of average temperatures of each study and save them

# Create a directory for output files
output_dir <- "data-processed/timetemps_df"

if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}


## Study 10: Alternating temperatures cycled every 24 hrs ----
dtrs <- (dtr %>% filter(study_id == "10"))$dtr

for (i in dtrs) {
	
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/study10_dtr", round(i,0), ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study 119: Sinusoidal fluctuations w/ a 3h period ----
dtr %>% filter(study_id == "119")

timetemps_df_119_4 <- timetemps_sin(dtr = 4, period = 3, time_resolution = 1) # interval: every hr

ggplot(timetemps_df_119_4, aes(x = time_hr, y = Tavg_26)) +
	geom_line() + scale_x_continuous(breaks = seq(0, 24, by = 3)) + theme_bw()

## Remove the time_hr column when saving the temperature profile
write_csv(as.data.frame(timetemps_df_119_4[,-1]), paste0(output_dir, "/study119_dtr4.csv"))


## Study 13: See Fig. 1; I have digitized the time series of 18-34ºC (mean 26ºC) and 18-28ºC (mean 22ºC) ----
# Load the digitized time series
timetemps_df_13_16_obs <- read_csv("data-processed/timetemps_df/study13_dtr16_obs.csv")
timetemps_df_13_10_obs <- read_csv("data-processed/timetemps_df/study13_dtr10_obs.csv")

ggplot(timetemps_df_13_16_obs, aes(x = time, y = temp)) +
	geom_line()

ggplot(timetemps_df_13_10_obs, aes(x = time, y = temp)) +
	geom_line()

# Use data between day 1 and 2 only (24 hrs)
timetemps_df_13_16_obs <- timetemps_df_13_16_obs %>% 
	filter(time >= 1 & time <= 2.03)

timetemps_df_13_10_obs <- timetemps_df_13_10_obs %>% 
	filter(time >= 1 & time <= 2.02)

## Convert time to days and change the first time point to 0
timetemps_df_13_16_obs$time <- (timetemps_df_13_16_obs$time - timetemps_df_13_16_obs$time[1]) * 24
timetemps_df_13_10_obs$time <- (timetemps_df_13_10_obs$time - timetemps_df_13_10_obs$time[1]) * 24

# Generate the data frame
timetemps_df_13_16 <- timetemps_timeseries(time_series_df = timetemps_df_13_16_obs, 
										   #mean_temp = 26, 
										   time_resolution = 1)

timetemps_df_13_10 <- timetemps_timeseries(time_series_df = timetemps_df_13_10_obs, 
										   #mean_temp = 22, 
										   time_resolution = 1)

# Check the results
ggplot() +
	geom_line(data = timetemps_df_13_16, aes(x = time_hr, y = Tavg_26)) +
	geom_point(data = timetemps_df_13_16_obs, aes(x = time, y = temp))

ggplot() +
	geom_line(data = timetemps_df_13_10, aes(x = time_hr, y = Tavg_24)) +
	geom_point(data = timetemps_df_13_10_obs, aes(x = time, y = temp))

## Remove the time_hr column when saving the temperature profile
write_csv(timetemps_df_13_16[,-1], paste0(output_dir, "/study13_dtr16.csv"))
write_csv(timetemps_df_13_10[,-1], paste0(output_dir, "/study13_dtr10.csv"))


# According to the time series, the average temperatures are: 
mean(timetemps_df_13_16_obs$temp) # 26.6 ºC
mean(timetemps_df_13_10_obs$temp) # 23.7 ºC
## Update the mean_temp for study 13 in data
data$mean_temp[data$study_id == "13" & data$dtr == "10"] <- 23.7
data$mean_temp[data$study_id == "13" & data$dtr == "16"] <- 26.6


## Study 15: Alternating temperatures cycled every 12 hrs ----
dtrs <- (dtr %>% filter(study_id == "15"))$dtr

for (i in dtrs) {
	
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/study15_dtr", i, ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study 159: Alternating temperatures with two daily 12-h cycles ----
dtrs <- (dtr %>% filter(study_id == "159"))$dtr

for (i in dtrs) {
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/study159_dtr", i, ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study 165: 8 hrs at 30ºC, 16 hrs at 22ºC ----
timetemps_df_165_8_long <- data.frame(time = seq(0, 24, by = 1), 
									  temp = c(rep(30, 8), 
									  		 rep(22, 16), 
									  		 30) # the time series should start and end at the same temperature
									  ) 

# Generate the data frame
timetemps_df_165_8  <- timetemps_timeseries(time_series_df = timetemps_df_165_8_long, 
											mean_temp = 25, 
											time_resolution = 1)

write_csv(timetemps_df_165_8[,-1], paste0(output_dir, "/study165_dtr8.csv"))



## Study 168: see fig. 1 ----
(dtr %>% filter(study_id == "168"))$dtr

## Fig. 1 only shows temp fluc regime of DTR = 8 and 10 (and 6 but where does that come from?)
# Load the digitized time series
timetemps_df_168_10_obs <- read_csv("data-processed/timetemps_df/study168_dtr10_obs.csv")
timetemps_df_168_8_obs <- read_csv("data-processed/timetemps_df/study168_dtr8_obs.csv")

head(timetemps_df_168_10_obs)

## Change the time so that it starts at time 0
timetemps_df_168_10_obs$time <- timetemps_df_168_10_obs$time - timetemps_df_168_10_obs$time[1]
timetemps_df_168_8_obs$time <- timetemps_df_168_8_obs$time - timetemps_df_168_8_obs$time[1]

timetemps_df_168_10 <- timetemps_timeseries(time_series_df = timetemps_df_168_10_obs, time_resolution = 1)
timetemps_df_168_8 <- timetemps_timeseries(time_series_df = timetemps_df_168_8_obs, time_resolution = 1)


# Check the results
ggplot() +
	geom_line(data = timetemps_df_168_10, aes(x = time_hr, y = Tavg_5), colour = "blue") +
	geom_point(data = timetemps_df_168_10_obs, aes(x = time, y = temp))

ggplot() +
	geom_line(data = timetemps_df_168_8, aes(x = time_hr, y = Tavg_4), colour = "blue") +
	geom_point(data = timetemps_df_168_8_obs, aes(x = time, y = temp))

## Remove the time_hr column when saving the temperature profile
write_csv(timetemps_df_168_10[,-1], paste0(output_dir, "/study168_dtr10.csv"))
write_csv(timetemps_df_168_8[,-1], paste0(output_dir, "/study168_dtr8.csv"))


## Study 9: Alternating temperatures cycled every 24 hrs ----
dtrs <- (dtr %>% filter(study_id == "9"))$dtr

for (i in dtrs) {
	
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/study9_dtr", i, ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study 90: Alternating temperatures cycled every 24 hrs ----
dtrs <- (dtr %>% filter(study_id == "90"))$dtr

for (i in dtrs) {
	
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/study90_dtr", i, ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study VB19-1: Sinusoidal fluctuations w/ a 24h period ----
dtr %>% filter(study_id == "VB19-1")

timetemps_df_VB19.1_10 <- timetemps_sin(dtr = 10, period = 24, time_resolution = 1) # interval: every hr

## check
ggplot(timetemps_df_VB19.1_10, aes(x = time_hr, y = Tavg_25)) +
	geom_line()

## Remove the time_hr column when saving the temperature profile
write_csv(as.data.frame(timetemps_df_VB19.1_10[,-1]), paste0(output_dir, "/studyVB19-1_dtr10.csv"))



## Study VB19-2: 16 hrs at higher temp, 8 hrs at lower temp ----
## 15.6-25.6ºC; DTR = 10
timetemps_df_VB19.2_10_obs <- data.frame(
	time = seq(0, 24, by = 1), 
	temp = c(rep(25.6, 16), rep(15.6, 8), 25.6) # the time series should start and end at the same temperature
) 

## 21.1-32.2ºC; DTR = 11.1
timetemps_df_VB19.2_11_obs <- data.frame(
		time = seq(0, 24, by = 1), 
		temp = c(rep(32.2, 16), rep(21.1, 8), 32.2) # the time series should start and end at the same temperature
	) 

## 10-23.9ºC; DTR = 13.9
timetemps_df_VB19.2_14_obs <- data.frame(
	time = seq(0, 24, by = 1), 
	temp = c(rep(23.9, 16), rep(10, 8), 23.9) # the time series should start and end at the same temperature
) 

# Generate the data frame
timetemps_df_VB19.2_10  <- timetemps_timeseries(time_series_df = timetemps_df_VB19.2_10_obs, 
											time_resolution = 1)

timetemps_df_VB19.2_11  <- timetemps_timeseries(time_series_df = timetemps_df_VB19.2_11_obs, 
												time_resolution = 1)

timetemps_df_VB19.2_14  <- timetemps_timeseries(time_series_df = timetemps_df_VB19.2_14_obs, 
												time_resolution = 1)

## Remove the time_hr column when saving the temperature profile
write_csv(timetemps_df_VB19.2_10[,-1], paste0(output_dir, "/studyVB19-2_dtr10.csv"))
write_csv(timetemps_df_VB19.2_11[,-1], paste0(output_dir, "/studyVB19-2_dtr11.csv"))
write_csv(timetemps_df_VB19.2_14[,-1], paste0(output_dir, "/studyVB19-2_dtr14.csv"))


## Study VB19-3: Alternating temperatures with two daily 12-h cycles ----
dtrs <- (dtr %>% filter(study_id == "VB19-3"))$dtr

for (i in dtrs) {
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/studyVB19-3_dtr", round(i, 0), ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}


## Study VB19-4: 16 hrs at 29.4°C, 8 hrs at 18.3°C ----

timetemps_df_VB19.4_11_obs <- data.frame(
	time = seq(0, 24, by = 1), 
	temp = c(rep(29.4, 16), rep(18.3, 8), 29.4) # the time series should start and end at the same temperature
) 

# Generate the data frame
timetemps_df_VB19.4_11  <- timetemps_timeseries(time_series_df = timetemps_df_VB19.4_11_obs,
												mean_temp = (16*29.4+8*18.3)/24,
												time_resolution = 1)

# Check the results
ggplot() +
	geom_line(data = timetemps_df_VB19.4_11, aes(x = time_hr, y = Tavg_26), colour = "blue") +
	geom_line(data = timetemps_df_VB19.4_11_obs, aes(x = time, y = temp))

## Remove the time_hr column when saving the temperature profile
write_csv(timetemps_df_VB19.4_11[,-1], paste0(output_dir, "/studyVB19-4_dtr11.csv"))


## Study VB19-5: Alternating temperatures with two daily 12-h cycles ----
dtrs <- (dtr %>% filter(study_id == "VB19-5"))$dtr

for (i in dtrs) {
	# Generate the temp profiles
	timetemps_df <- timetemps_alt(dtr = i)
	
	# Define file path
	file_name <- paste0(output_dir, "/studyVB19-5_dtr", round(i, 0), ".csv")
	
	# Save the plot
	write_csv(as.data.frame(timetemps_df), file_name)
	print(paste("Saved:", file_name))
}



##########
###### 2. Rate summation calculations ----
##########

## Create a data frame that show the study_id of each curves and the DTR tested in that study_id
info <- data %>% 
	filter(dtr != 0) %>% 
	group_by(curve_id, dtr) %>% 
	distinct(curve_id, study_id, dtr)

## Round the dtr
info$dtr <- round(info$dtr, 0)
head(info)

## Since I cannot find the time series for DTR = 4 and DTR = 9 for study168, remove those
info <- info %>% 
	filter(study_id != "168" | 
		   	study_id == "168" & !(dtr %in% c(4,9)))


## Get all the curve_ids
curve_ids <- unique(info$curve_id)
curve_ids

# Create a directory to save RS calculations
output_dir <- "data-processed/RS-calculations"

if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}


# Add in a progress bar because the following code takes ~40 mins
pb <- progress_bar$new(
	format = "[:bar] :percent eta: :eta",
	total = nrow(info)
)

## loop over each curve_id
for (i in curve_ids){
	# Get the study_id of that curve_id
	study_id <- pull(d_preds_best %>% filter(curve_id == i) %>% distinct(study_id))
	
	# Get a list of DTR's tested in study_id
	dtrs <- pull(info %>% filter(curve_id == i), 3)
	
	# Loop over each DTR tested in the study
	for (k in dtrs) {
		message(paste0("Fitting: curve_id = ", i,
					 ", study_id = ", study_id,
					 ", DTR = ", k))
		
		# Increment progress bar
		pb$tick()
		
		## Wrap in tryCatch to skip errors
		tryCatch({
			## Load the temperature time series of that study_id and DTR
			timetemps_df <- read.csv(paste0("data-processed/timetemps_df/study", study_id, "_dtr", k, ".csv"), 
									 check.names = FALSE, stringsAsFactors = FALSE)
			timetemps_df <- as.data.frame(timetemps_df)
			
			## Get the TPC model fit of that curve_id
			model_fit <- d_preds_best %>% filter(curve_id == i)
			
			## Apply rate summation 
			pred_RS <- RSCalcTempGrad(model_fit, timetemps_df, temp_grad)
			
			## Save the results as a csv file
			file_name <- paste0("data-processed/RS-calculations/curve", i, "_dtr", k, ".csv")
			write_csv(pred_RS, file_name)
		},
		error = function(e) {
			message(paste0("⚠️ Error for curve_id = ", i,
						   ", study_id = ", study_id, 
						   ", DTR = ", k, 
						   ": ", e$message))
		})
	}
}

## debug section ----
## Curves 4-8, 30, 39 (DTR8, 10), 40(8,10), 42(8,10), 48(10), 49(5,10), 53(12, 10), 65,66,67

# dtr_test <- dtr %>% subset(study_id == "9" | study_id == "10")
# curve_ids <- unique(test$curve_id)
# curve_ids

test <- d_preds_best %>% filter(curve_id == 4)

timetemps_df <- read.csv(paste0("data-processed/timetemps_df/study13_dtr10.csv"), 
						 check.names = FALSE, stringsAsFactors = FALSE)

TPC_predictions <- d_preds_best %>% filter(curve_id == 4)
pred_RS <- RSCalcTempGrad(TPC_predictions, timetemps_df, temp_grad) # error, let's break down the function

predictions_long <- TPC_predictions %>%
	select(temp, .fitted, conf_lower, conf_upper) %>% 
	pivot_longer(!temp, names_to = "type", values_to = "trait_value")  %>%
	mutate(temp = as.numeric(temp)) %>% # make temperature numeric
	relocate(type)

rs_calc_gradient <- as.data.frame(matrix(nrow = 3, #.fitted, conf_lower, conf_upper
										 ncol = ncol(timetemps_df)))

i = 290
temp_seq <- data.frame(temp = timetemps_df[, i])

sj <- left_join(temp_seq, predictions_long, by = "temp", relationship = "many-to-many")

if (any(is.na(sj$trait_value))) {
	rs_calc_gradient[, i] <- NA
} else {
	rs_calc <- sj %>%
		group_by(type) %>%
		summarise(rs_calc = mean(trait_value))
	
	# Store rate summation calculation in the corresponding column
	rs_calc_gradient[, i] <- pull(rs_calc[, 2])
}


  rs_calc_gradient <- as.data.frame(t(rs_calc_gradient))
  
  colnames(rs_calc_gradient) <- c("mean", "conf_lower", "conf_upper")
  
  rs_calc_gradient <- rs_calc_gradient %>% 
  	mutate(temp = temp_grad) %>% 
  	relocate(temp)

### debug section end ----

##########
###### 3. Plotting ----
##########

## Get all the curve_ids
curve_ids <- unique(info$curve_id)
curve_ids

## # Create a directory to save the plots
output_dir <- "figures/Lilian/tpc-fitted-20250618/rs"

if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}

# Define custom colors for each model
model_colors <- c(
	"const_tpc" = "#000075",
	"rs" = "#EFC000FF",
	"const" = "black",
	"fluc" = "firebrick1")

# Define custom labels
model_labels <- c(
	"const_tpc" = "Constant TPC",
	"rs" = "RS predictions",
	"const" = "Constant",
	"fluc" = "Fluctuation")

	
# Add in a progress bar because the following code
pb <- progress_bar$new(
	format = "[:bar] :percent eta: :eta",
	total = nrow(info)
)

## loop over each curve_id
for (i in curve_ids){
	# Get the study_id of that curve_id
	study_id <- pull(info %>% filter(curve_id == i), 2)
	
	# Get a list of DTR's tested in study_id
	dtrs <- pull(info %>% filter(curve_id == i), 3)
	
	# Loop over each DTR tested in the study
	for (k in dtrs) {
		# Increment progress bar
		pb$tick()
		
		# Create subset of data
		## raw data
		### The round ensures to even if dtr is not exactly as the one in info we can still get the raw data
		mini_data <- data %>% filter(curve_id == i) %>% filter(dtr == 0 | round(dtr,0) == k)
		## Constant temp TPC
		mini_pred <- d_preds_best %>% filter(curve_id == i)
		## RS calculation
		mini_rs <- read.csv(paste0("data-processed/RS-calculations/curve", i, "_dtr", k, ".csv"), 
							check.names = FALSE, stringsAsFactors = FALSE)
		
		# Create a ggplot
		RSplot <- ggplot() +
			# Constant temp TPC
			geom_ribbon(data = mini_pred, 
						aes(x = temp, ymin = conf_lower, ymax = conf_upper), 
						fill = "#4363d8", 
						alpha = 0.5
						) +
			geom_line(data = mini_pred, 
					  aes(x = temp, y = .fitted, 
					  	color = "const_tpc"), 
					  linewidth = 1
					  ) +
			
			# Rate summation predictions
			geom_ribbon(data = mini_rs,
						aes(x = temp, ymin = conf_lower, ymax = conf_upper), 
						fill = "gold1", 
						alpha = 0.5
						) +
			geom_line(data = mini_rs,
					  aes(x = temp, y = mean, color = "rs"), 
					  linewidth = 1
					  ) +
			
			# Raw data
			geom_point(data = mini_data, 
					   aes(x = mean_temp, y = response, 
					   	color = ifelse(dtr == 0, "const", "fluc")), 
					   size = 2
					   ) +
			geom_hline(yintercept = 0, linetype = "dashed") +
			ylim(0, NA) +
			
			# Customize labels
			labs(title = paste("CurveID", i, 
							   "StudyID", unique(mini_data$study_id), 
							   "DTR", k, unique(mini_pred$model_name), 
							   unique(mini_data$trait), 
							   unique(mini_data$additional_complexity)),
				 x = expression(paste("Temperature (", degree, "C)")),
				 y = paste(unique(mini_data$trait), "(", mini_data$units, ")")) +
			 scale_color_manual(name = "",
			 				   values = model_colors, labels = model_labels,
			 				   breaks = c("const", "fluc", "const_tpc", "rs")
			 				   ) +
			theme_bw()
		
		# Define file path
		file_name <- paste0(output_dir, "/plot_", i, "_", 
							unique(mini_pred$study_id), "_", k, "_rs.png")
		
		# Save the plot
		ggsave(file_name, plot = RSplot, width = 8, height = 4, dpi = 300)
		
		message(paste("Saved:", file_name))
	}
}

