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
#### Function to perform rate summation ----
#### Written by Marta Shocket in 2022; See: https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/9f15762d8a083f45a58656a15a9c5eb49bf02da6/R-scripts/working-versions-code/00_RSProjectFunctions.R

# Arguments:
#	1) TPC_predictions = predicted trait values (from constant temperature TPC) across the temperature gradient for every MCMC iteration 
#			df: 451 cols (temperature gradient) x 5000 rows (MCMC iterations)
# 2) timetemps_df = hourly temperatures across the mean temperature gradient for a given dtr treatment
#			df: 451 cols (mean temperature gradient) x 24 rows (24 hourly temperature values)
# 3) temp_grad = temperature gradient values to apply rate summation over
#			vector: 451 elements
#
# Returns: 
#	1) rs_calc_gradient = trait values predicted by rate summation across a mean temperature gradient for every MCMC iteration
#			df: 451 cols (temperature gradient) x 5000 rows (MCMC iterations)
#

RSCalcTempGrad <- function(TPC_predictions, timetemps_df, temp_grad) {
  
  # Name the columns so the pivot_longer() below will work
  colnames(TPC_predictions) <- temp_grad
  
  # Pivot constant temperature TPC trait predictions into long format:
  #		3 cols (iteration, temperature, & trait_value) x 3.4 million rows (7500 MCMC iterations x 451 temperatures)
  # Note: this is the longest step, do it as few times as possible
  predictions_long <- TPC_predictions %>%
    mutate(iteration = rownames(.)) %>% # add column with iteration number 
    pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value")  %>%
    mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
    mutate(iteration = as.numeric(iteration))
  
  # Create empty data frame to store rate summation calculations
  rs_calc_gradient <- as.data.frame(matrix(nrow = nrow(TPC_predictions), ncol = ncol(timetemps_df)))
  
  # Add in a progress bar because this function takes a way
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = ncol(timetemps_df)
  )
  
  # Loop through each mean temperature along the gradient
  for (i in 1:ncol(timetemps_df)) {
    
    # Increment progress bar
    pb$tick()
    
    # Pull out the hourly temperatures for that mean temperature as a df with one col
    # NOTE: the column name must match the pivoted trait predictions column for the left_join to work
    temp_seq <- data.frame(temperature = timetemps_df[, i])
    
    # Left join temperature sequence and trait prediction:
    #		for every row in temp_seq, it adds every trait value + iteration combo from predictions_long with that temperature value
    #		-> df with 3 cols (temperature, iteration, trait value) and 180,000 rows (7500 MCMC x 24 hourly temperatures)
    sj <- left_join(temp_seq, predictions_long, by = "temperature", relationship = "many-to-many")
    
    # Rate summation calculation averages trait values at all hourly temperatures for a given MCMC iteration
    #		-> df with 2 cols (iteration, rs_calc) and 7500 rows (MCMC)
    rs_calc <- sj %>%
      group_by(iteration) %>%
      summarise(rs_calc = mean(trait_value))
    
    # Store rate summation calculation in the corresponding column
    rs_calc_gradient[, i] <- pull(rs_calc[, 2])
    
  } # End loop through mean temperature gradient
  
  return(rs_calc_gradient)
  
} # End function

