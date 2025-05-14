# Define a function to calculate trait value under fluctuating temperature
rate_summation <- function(mean_temp, tpc){
  temp_max <- mean_temp + 5
  temp_min <- mean_temp - 5
  
  # Extract the trait values from the TPC 
  trait_max <- tpc %>%
    filter(temp == temp_max) %>%
    pull(mean)
  
  trait_min <- tpc %>% 
    filter(temp == temp_min) %>% 
    pull(mean)
  
  # Debugging outputs
  cat("Mean Temp:", mean_temp, "\n")
  cat("Temp Max:", temp_max, "Trait Max:", trait_max, "\n")
  cat("Temp Min:", temp_min, "Trait Min:", trait_min, "\n")
  
  # Calculate the average trait values for fluctuating temperature
  # If missing one value, then return NA
  fluct_trait <- mean(c(trait_max, trait_min), na.rm = F) 
  cat("Trait value under fluctuating temp:", fluct_trait, "\n\n")
  return(fluct_trait)
}

rate_summation_min <- function(mean_temp, tpc){
  temp_max <- mean_temp + 5
  temp_min <- mean_temp - 5
  
  # Extract the trait values from the TPC 
  trait_max_ci <- tpc %>%
    filter(temp == temp_max) %>%
    pull(X2.5.)
  
  trait_min_ci <- tpc %>% 
    filter(temp == temp_min) %>% 
    pull(X2.5.)
  
  # Debugging outputs
  cat("Mean Temp:", mean_temp, "\n")
  cat("Temp Max:", temp_max, "Trait Max:", trait_max_ci, "\n")
  cat("Temp Min:", temp_min, "Trait Min:", trait_min_ci, "\n")
  
  # Calculate the average trait values for fluctuating temperature
  # If missing one value, then return NA
  fluct_trait <- mean(c(trait_max_ci, trait_min_ci), na.rm = F) 
  cat("Trait value under fluctuating temp:", fluct_trait, "\n\n")
  return(fluct_trait)
}

rate_summation_max <- function(mean_temp, tpc){
  temp_max <- mean_temp + 5
  temp_min <- mean_temp - 5
  
  # Extract the trait values from the TPC 
  trait_max_ci <- tpc %>%
    filter(temp == temp_max) %>%
    pull(X97.5.)
  
  trait_min_ci <- tpc %>% 
    filter(temp == temp_min) %>% 
    pull(X97.5.)
  
  # Debugging outputs
  cat("Mean Temp:", mean_temp, "\n")
  cat("Temp Max:", temp_max, "Trait Max:", trait_max_ci, "\n")
  cat("Temp Min:", temp_min, "Trait Min:", trait_min_ci, "\n")
  
  # Calculate the average trait values for fluctuating temperature
  # If missing one value, then return NA
  fluct_trait <- mean(c(trait_max_ci, trait_min_ci), na.rm = F) 
  cat("Trait value under fluctuating temp:", fluct_trait, "\n\n")
  return(fluct_trait)
}