## Rate Summation Project
## Functions to load for use in other scripts

## Table of contents:
## 1. Function to perform rate summation


#### 0. Load tidyverse library ----
library(tidyverse)

################################################################################
#### 1. Functions to process TPC model output
######## A. Function to calculate TPC posterior summary statistics across a temperature gradient

calcPostQuants <- function(TPC_predictions, trait_treatment_name, temp_gradient) {
  
  # Reassign column names to the temperature gradient
  colnames(TPC_predictions) <- temp_gradient
  
  output <- TPC_predictions %>% 
    mutate(iteration = rownames(.)) %>% # add column with iteration number 
    pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
    group_by(temperature) %>%
    summarise(lowerCI = quantile(trait_value, probs = 0.025),
              upperCI = quantile(trait_value, probs = 0.975),
              mean = mean(trait_value),
              median = median(trait_value)) %>% 
    mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
    arrange(temperature) %>% # re-order rows by ascending temperature (since grouping made it categorical, it is ordered alphabetical)
    mutate(treatment = trait_treatment_name) # add column with variable + treatment name
  
  return(output) # return output
  
}


################################################################################
#### Function to do Bootstrapping ----

# Arguments:
# 1) TPC_fits = predicted trait values (from constant temperature TPC) across the temperature gradient for every MCMC iteration 
#			df: 451 cols (temperature gradient) x 5000 rows (MCMC iterations)
# 2) model_name = the TPC model fitted. Must be "briere", "quadratic", "lactin", or "sharpeschoolfull"
#			
# 3) temp_grad = temperature gradient values to apply rate summation over
#			vector: 451 elements
#
# Returns: 
#	1) rs_calc_gradient = trait values predicted by rate summation across a mean temperature gradient for every MCMC iteration
#			df: 451 cols (temperature gradient) x 5000 rows (MCMC iterations)
#

bootstrapping <- function(TPC_fits, model_name, temp_grad) {
	# run for loops to bootstrap the TPC model
	for(i in 1:nrow(TPC_fits)){
		temp_data <- TPC_fits$data[[i]]
		temp_fit <- nlsLM(rate ~ briere1_1999(temp = temp, tmin, tmax, a),
						  data = temp_data,
						  start = TPC_fits$coefs[[i]],
						  lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'),
						  upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'))
		boot <- Boot(temp_fit, method = 'residual')
		TPC_fits$bootstrap[[i]] <- boot
		rm(list = c('temp_fit', 'temp_data', 'boot'))
	}
	
	## Get the raw values of each bootstrap
	TPC_fits <- mutate(TPC_fits, output_boot = map(bootstrap, function(x) x$t))
	
	## Calculate predictions with a gnarly written function
	TPC_fits <- mutate(TPC_fits, preds = map2(output_boot, data, function(x, y){
		temp <- as.data.frame(x) %>%
			drop_na() %>%
			mutate(iter = 1:n()) %>%
			group_by_all() %>%
			do(data.frame(temp = temp_grad)) %>%
			ungroup() %>%
			mutate(pred = briere1_1999(temp = temp, tmin, tmax, a))
		return(temp)
	}))
	
	# select, unnest and calculate 95% CIs of predictions
	boot_conf_preds_bri <- select(TPC_fits, curve_id, preds) %>%
		unnest(preds) %>%
		drop_na(pred) %>% 
		group_by(curve_id, temp) %>%
		summarise(conf_lower = quantile(pred, 0.025),
				  conf_upper = quantile(pred, 0.975),
				  .groups = 'drop') %>% 
		mutate(model_name = model_name)
}





################################################################################
#### Function to perform rate summation ----
#### The following code was adapted from https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/9f15762d8a083f45a58656a15a9c5eb49bf02da6/R-scripts/working-versions-code/00_RSProjectFunctions.R

# Arguments:
# 1) TPC_predictions = predicted trait values (from constant temperature TPC) across the temperature gradient for every MCMC iteration 
#			df: 501 columns (temperature gradient from 0 to 50ºC at 0.1ºC interval) x 3 rows (.fitted, conf_lower, conf_upper)
# 2) timetemps_df = hourly temperatures across the mean temperature gradient for a given dtr treatment
#			df: 501 cols (mean temperature gradient) x 24 rows (24 hourly temperature values)
# 3) temp_grad = temperature gradient values to apply rate summation over
#			vector: 501 elements
#
# Returns: 
# 1) rs_calc_gradient = trait values predicted by rate summation across a mean temperature gradient for every MCMC iteration
#			df: 501 cols (temperature gradient) x 5000 rows (MCMC iterations)
#

RSCalcTempGrad <- function(TPC_predictions, timetemps_df, temp_grad) {
  
  # Pivot constant temperature TPC trait predictions into long format:
  #		3 cols (type, temp, & trait_value) x 3.4 million rows (7500 MCMC iterations x 451 temperatures)
  # Note: this is the longest step, do it as few times as possible
  predictions_long <- TPC_predictions %>%
  	select(temp, .fitted, conf_lower, conf_upper) %>% 
  	pivot_longer(!temp, names_to = "type", values_to = "trait_value")  %>%
  	mutate(temp = as.numeric(temp)) %>% # make temperature numeric
  	relocate(type)
  
  # Create empty data frame to store rate summation calculations
  rs_calc_gradient <- as.data.frame(matrix(nrow = 3, #.fitted, conf_lower, conf_upper
  										   ncol = ncol(timetemps_df)))
  
  # Add in a progress bar because this function takes a way
  # pb <- progress_bar$new(
  #   format = "[:bar] :percent eta: :eta",
  #   total = ncol(timetemps_df)
  # )
  
  # Loop through each mean temperature along the gradient
  for (i in 1:ncol(timetemps_df)) {
    
    # Increment progress bar
    #pb$tick()
    
    # Pull out the temperature time series for that mean temperature as a df with one col
    # NOTE: the column name must match the pivoted trait predictions column for the left_join to work
    temp_seq <- data.frame(temp = timetemps_df[, i])
    
    # Left join temperature sequence and trait prediction:
    #		for every row in temp_seq, it adds every trait value + iteration combo from predictions_long with that temperature value
    #		-> df with 3 cols (temperature, iteration, trait value) and 180,000 rows (7500 MCMC x 24 hourly temperatures)
    sj <- left_join(temp_seq, predictions_long, by = "temp", relationship = "many-to-many")
    
    # Rate summation calculation averages trait values at all hourly temperatures for a given MCMC iteration
    #		-> df with 2 cols (type, rs_calc) and 3 rows (.fitted, conf_lower, conf_upper)
    # If there are any NAs in trait_value, fill entire column with NA
    if (any(is.na(sj$trait_value))) {
    	rs_calc_gradient[, i] <- NA
    } else {
    	rs_calc <- sj %>%
    		group_by(type) %>%
    		summarise(rs_calc = mean(trait_value))
    	
    	# Store rate summation calculation in the corresponding column
    	rs_calc_gradient[, i] <- pull(rs_calc[, 2])
    }
    
  } # End loop through mean temperature gradient
  
  ## Change the rs_calc_gradient data frame into a different format
  rs_calc_gradient <- as.data.frame(t(rs_calc_gradient))
  
  colnames(rs_calc_gradient) <- c("mean", "conf_lower", "conf_upper")
  
  rs_calc_gradient <- rs_calc_gradient %>% 
  	mutate(temp = temp_grad) %>% 
  	relocate(temp)
  
} # End function



# Function to generate a data frame of sinusoidal temperature fluctuations ----
# across a range of average temperatures and a full 24-hour period

## Arguments:
## 1) temp_grad = sequence of average temperatures from 0 to 50ºC at 0.1ºC interval
##			vector: 501 elements
## 2) dtr = diurnal temperature range (maximum temperature - minimum temp)
## 3) period = # Period of sinusoidal cycle in hours
## 4) time_resolution = time step interval (in hours)
##
##
## Returns: 
##	1) timetemps_df = data frame with the temperature profiles over time in a sinusoidal temperature fluctuating regime, for a range of average temperatures
##			df: 501 cols (average temperatures from 0°C to 50°C in 0.1°C 
##              increments) x 24/time_resolution rows (Time steps over 24 hours)
##

timetemps_sin <- function(
		temp_grad = seq(0, 50, 0.1), # Sequence of average temperatures from 0ºC to 50ºC
		dtr,                         # 2*Amplitude of the temperature fluctuation (°C)
		period,                      # Period of sinusoidal cycle in hours
		time_resolution = 1       # Time resolution in hours (Default: 1 hr)
) {
	
	# Generate time steps from 0 to 24 hours with specified resolution
	time_steps <- seq(0, 23, by = time_resolution)
	
	# Initialize an empty matrix to hold temperature values
	temp_matrix <- matrix(nrow = length(time_steps), ncol = length(temp_grad))
	
	# Loop over each average temperature and calculate sinusoidal values
	for (i in seq_along(temp_grad)){
		T_avg <- temp_grad[i]
		
		# Sinusoidal temperature: T(t) = T_avg + A * sin(2π / P * t), and A = dtr/2
		temp_matrix[, i] <- T_avg + (dtr/2) *sin(2*pi / period * time_steps)
	}
	
	# Convert the matrix to a data frame
	timetemps_df <- as.data.frame(temp_matrix)
	
	# Assign the average temperatures to be column names
	colnames(timetemps_df) <- paste0("Tavg_", temp_grad)
	
	# Add a column to show the time step
	timetemps_df$time_hr <- time_steps
	timetemps_df <- timetemps_df %>% relocate(time_hr) # move time column to the front
	
	# Change the temperature value to 0 if it is below 0 (to ensure they are within the model prediction)
	timetemps_df[timetemps_df < 0] <- 0
	# Change the temperature value to 0 if it is above 50
	timetemps_df[timetemps_df > 50] <- 50
	
	# Round to nearest 0.1C and convert to df
	timetemps_df = round(as.data.frame(timetemps_df), 1)
	return(timetemps_df)
}


# Function to generate a data frame of alternating temperature fluctuations ----
# across a range of average temperatures and a full 24-hour period

## Arguments:
## 1) temp_grad = sequence of average temperatures from 0 to 50ºC at 0.1ºC interval
##			vector: 501 elements
## 2) dtr = diurnal temperature range (maximum temperature - minimum temp)
##
##
## Returns: 
##	1) timetemps_df = data frame with the temperature profiles in an alternating fluctuation regime over time, for a range of average temperatures
##			df: 501 cols (average temperatures from 0°C to 50°C in 0.1°C 
##              increments) x 2 rows (min and max temp)
##

timetemps_alt <- function(temp_grad = seq(0, 50, 0.1), dtr) {
	timetemps_df <- data.frame(min = round(temp_grad - dtr/2, 1), max = round(temp_grad + dtr/2, 1))
	
	# Change into wide format
	timetemps_df <- as.data.frame(t(timetemps_df))
	
	# Change the temperature value to 0 if it is below 0 (to ensure they are within the model prediction)
	timetemps_df[timetemps_df < 0] <- 0
	# Change the temperature value to 0 if it is above 50
	timetemps_df[timetemps_df > 50] <- 50
	
	# Assign the average temperatures to be column names
	colnames(timetemps_df) <- paste0("Tavg_", temp_grad)
	
	# Round to nearest 0.1C and convert to df
	timetemps_df = round(as.data.frame(timetemps_df), 1)
	
	return(timetemps_df)
}


# Function to generate a data frame of other diurnal temperature fluctuation across a range of average temperatures ----

## This function takes a digitized time series (e.g., extracted from a figure)
## representing the thermal fluctuations used in the experiment over 24 hrs and 
## scales that fluctuation pattern across a range of average temperatures
##
## Arguments:
## 1) temp_grad = sequence of average temperatures from 0 to 50ºC at 0.1ºC interval
##			vector: 501 elements
## 2) time_series_df = temperature time series around one known average temperature, start and end at the same temperature
##			df: 2 cols (time: time point in hrs, temp: temperature at the corresponding time point)
## 3) mean_temp = average temperature of the time series
##			if mean_temp is not provided, it will be the average temp over the time series
## time_resolution = time resolution in hours (Default: 1 hr; if time resolution is 30 mins, then time_resolution = 0.5)
##
##
## Returns: 
##	1) timetemps_df = data frame with the temperature profiles in an diurnal fluctuation regime over time, for a range of average temperatures
##			df: 501 cols (average temperatures from 0°C to 50°C in 0.1°C 
##              increments) x 24/time_resolution rows (Time steps over 24 hours)

timetemps_timeseries <- function(
		temp_grad = seq(0, 50, 0.1), 
		time_series_df,              
		mean_temp = mean(time_series_df$temp),                   
		time_resolution = 1
		) {
	# Interpolate to evenly spaced time points
	approx_temp <- approx(x = time_series_df$time, y = time_series_df$temp, 
						  xout = seq(min(time_series_df$time), 
						  		   max(time_series_df$time), 
						  		   length.out = 24 / time_resolution + 1))$y
	
	# Center the interpolated time series around mean temperature
	centered_temp <- approx_temp - mean_temp
	# Generate a matrix: for each mean temp, shift the centered profile
	temp_matrix <- sapply(temp_grad, function(mu) centered_temp + mu)
	
	# Convert the matrix to a data frame
	timetemps_df <- as.data.frame(temp_matrix)
	
	# Assign the average temperatures to be column names
	colnames(timetemps_df) <- paste0("Tavg_", temp_grad)
	
	# Add a column to show the time step
	timetemps_df$time_hr <- seq(0, 24, by = time_resolution)
	timetemps_df <- timetemps_df %>% relocate(time_hr) # move time column to the front
	
	# Change the temperature value to 0 if it is below 0 (to ensure they are within the model prediction)
	timetemps_df[timetemps_df < 0] <- 0
	# Change the temperature value to 0 if it is above 50
	timetemps_df[timetemps_df > 50] <- 50
	
	
	# Remove the last row since time = 24 is equal to time = 0
	timetemps_df <- timetemps_df[-nrow(timetemps_df), ]
	
	# Round to nearest 0.1C and convert to df
	timetemps_df = round(as.data.frame(timetemps_df), 1)
	
	return(timetemps_df)
}
