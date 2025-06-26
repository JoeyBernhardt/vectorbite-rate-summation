## Lilian Chan, University of Guelph
## Rate summation projection
##
## Purpose: Fit TPCs (Briere, quadratic, Sharpe-Schoolfield (full), Lactin) 
## to the data from constant temperature treatment
## 
##
## Table of content:
##    0. Set-up workspace
##    1. Preparing the data
##    2. Rate summation calculations
##    3. Plotting


##########
###### 0. Set-up workspace ----
##########

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
###### 1. Preparing the data  ----
##########

## Load data
data <- read_csv("data-processed/vb_maggie_data_20250611.csv")

## Load the best-fitting constant temperature TPC of each curve
load("R-scripts/model-output/d_preds_best.Rdata")

## Load all constant temperature TPCs of each curve
load("R-scripts/model-output/d_preds.Rdata")

# Round the temperature to the nearest 0.1 (to avoid errors when doing joining later)
d_preds_best$temp <- round(d_preds_best$temp, 1)
d_preds$temp <- round(d_preds$temp, 1)

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


## Create a data frame that show the study_id of each curves and the DTR tested in that study_id
info <- data %>% 
	filter(dtr != 0) %>% 
	group_by(curve_id, dtr) %>% 
	distinct(curve_id, study_id, dtr)

info$dtr <- round(info$dtr, 0)


## Since I cannot find the time series for DTR = 4 and DTR = 9 for study168, remove those
info <- info %>% 
	filter(study_id != "168" | 
		   	study_id == "168" & !(dtr %in% c(4,9)))

head(info)
##########
###### 2. Rate summation calculations ----
##########

info <- info %>% filter(study_id == "150")

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

## loop over each curve_id and dtr in info
for (i in 1:nrow(info)){
	
	# Get the curve_id
	curveid <- pull(info[i, 1])
	
	# Get the study_id of that curve_id
	studyid <- pull(info[i, 2])
	
	# Get the DTR
	dtr <- pull(info[i, 3])
	
	message(paste0("Fitting: curve_id = ", curveid,
				   ", study_id = ", studyid,
				   ", DTR = ", dtr))
	
	# Increment progress bar
	pb$tick()
	
	## Wrap in tryCatch to skip errors
	tryCatch({
		## Load the temperature time series of that study_id and DTR
		timetemps_df <- read.csv(paste0("data-processed/timetemps_df/study", studyid, "_dtr", dtr, ".csv"), 
								 check.names = FALSE, stringsAsFactors = FALSE)
		timetemps_df <- as.data.frame(timetemps_df)
		
		## Get the TPC model fit of that curve_id
		model_fit <- d_preds_best %>% filter(curve_id == curveid)
		
		## Apply rate summation 
		pred_RS <- RSCalcTempGrad(model_fit, timetemps_df, temp_grad)
		
		## Save the results as a csv file
		file_name <- paste0("data-processed/RS-calculations/curve", curveid, "_dtr", dtr, ".csv")
		write_csv(pred_RS, file_name)
	},
	error = function(e) {
		message(paste0("⚠️ Error for curve_id = ",curveid ,
					   ", study_id = ", studyid, 
					   ", DTR = ", dtr, 
					   ": ", e$message))
	})
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
for (i in 1:nrow(info)){
	
	# Get the curve_id
	curveid <- pull(info[i, 1])
	
	# Get the study_id of that curve_id
	studyid <- pull(info[i, 2])
	
	# Get the DTR
	DTR <- pull(info[i, 3])
	
	# Increment progress bar
	pb$tick()
	
	# Create subset of data
	## raw data
	### The round ensures to even if dtr is not exactly as the one in info we can still get the raw data
	mini_data <- data %>% filter(curve_id == curveid) %>% filter(dtr == 0 | round(dtr,0) == DTR)
	## Constant temp TPC
	mini_pred <- d_preds_best %>% filter(curve_id == curveid)
	## RS calculation
	mini_rs <- read.csv(paste0("data-processed/RS-calculations/curve", curveid, "_dtr", DTR, ".csv"), 
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
		labs(title = paste("CurveID", curveid, 
						   "StudyID", unique(mini_data$study_id), 
						   "DTR", DTR, unique(mini_pred$model_name), 
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
	file_name <- paste0(output_dir, "/plot_", curveid, "_", 
						unique(mini_pred$study_id), "_", DTR, "_rs.png")
	
	# Save the plot
	ggsave(file_name, plot = RSplot, width = 8, height = 4, dpi = 300)
	
	message(paste("Saved:", file_name))	
	
}

