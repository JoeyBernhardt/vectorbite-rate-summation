## Lilian Chan, University of Guelph
## Rate summation projection
##
## Purpose: Fit TPCs (Briere, quadratic, Sharpe-Schoolfield (full), Lactin) 
## to the data from constant temperature treatment
## 
##
## Table of content:
##    0. Set-up workspace
##
##    1. Prepare the data
##    2. TPCs fitting using rTPC
##        A. Fitting
##           i.  Briere and Quadratic (3 parameters)
##          ii.  Lactin (4 parameters)
##         iii.  Sharpe-Schoolfield (6 parameters)
##        B. Predicting
##        C. Bootstrapping
##        D. Plotting



##########
###### 0. Set-up workspace ----
##########
## load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

setwd("~/Documents/UofG/vectorbite-rate-summation")

# load in data
data <- read_csv("data-processed/vb_maggie_data_20250507.csv")

data$dtr <- as.factor(data$dtr)

## Load the Quadratic function
source("R-scripts/00_Quadratic.R")


##########
###### 1. Prepare the data ----
##########

## Only fit TPC to the constant treatments data
data_const <- data %>% filter(dtr == 0)

## Change column names for mean_temp and response to prevent mismatch  
colnames(data_const)[6] <- "temp"
colnames(data_const)[18] <- "rate"

## Select columns we need
data_const <- data_const %>% 
  filter(!is.na(rate)) %>% 
  select(curve_id, study_id, trait, trait_def, units, genus, species, 
         additional_complexity, temp, rate)



##########
###### 2. TPCs fitting using rTPC ----
######     Ai. Fitting: Briere and Quadratic (3 parameters)  ----
##########

## The following code was adapted from https://padpadpadpad.github.io/rTPC/articles/fit_many_curves.html

## Edit nls_multstart to allow for a progress bar
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower, 
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
                                   control, modelweights, ...){
  if(!is.null(pb)){
    pb$tick()
  }
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower, 
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count, 
                control = control, modelweights = modelweights, ...)
}

## start progress bar and estimate time it will take
number_of_models <- 2
number_of_curves <- length(unique(data_const$curve_id))

## setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")


## fit the chosen model formulation in rTPC
d_fits <- nest(data_const, data = c(temp, rate)) %>%
  mutate(briere = map(data, 
                      ~nls_multstart_progress(
                        rate~briere1_1999(temp = temp, tmin, tmax, a),
                        data = .x,
                        iter = c(3,3,3),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere1_1999') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere1_1999') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere1_1999'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere1_1999'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)
                      ), 
         quadratic = map(data, 
                      ~nls_multstart_progress(
                        formula = rlang::expr(rate ~ quadratic(temp, tmin, tmax, a)), #rlang::expr is very very important. Without it the code below won't run
                        data = .x,
                        iter = c(3,3,3),
                        start_lower = quadratic.starting_vals(.x$temp, .x$rate) - 10,
                        start_upper = quadratic.starting_vals(.x$temp, .x$rate) + 10,
                        lower = quadratic.lower_lims(.x$temp, .x$rate),
                        upper = quadratic.upper_lims(.x$temp, .x$rate),
                        supp_errors = 'Y',
                        convergence_count = FALSE)
                      )
         )



##########
###### 2Aii. Fitting: Lactin (4 parameters)  ----
##########

## Since the Lactin model has 4 parameters, we need curves with >= 4 data points from constant temp treatment
num_of_const <- data %>% 
  filter(dtr == 0) %>% 
  group_by(curve_id) %>% 
  summarise(n = n()) %>% 
  arrange(n)

curve_needed_lac <- num_of_const %>% 
  filter(n >= 4)

curve_needed_lac <- c(curve_needed_lac$curve_id) 
curve_needed_lac # 60 curves with >= 4 const temp treatment


data_const_lac <- data_const %>% 
  filter(curve_id %in% curve_needed_lac)

unique(data_const_lac$curve_id) # should be 60 curves


## start progress bar and estimate time it will take
number_of_models <- 1
number_of_curves <- length(unique(data_const_lac$curve_id))

## setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")


## fit the chosen model formulation in rTPC
d_fits_lac <- nest(data_const_lac, data = c(temp, rate)) %>%
  mutate(lactin = map(data, 
                      ~nls_multstart_progress(
                        rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                        data = .x,
                        iter = c(3,3,3,3),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)
                      )
         )


#Combine d_fits and d_fits
colnames(d_fits)
d_fits <- full_join(d_fits, d_fits_lac[,c("curve_id", "lactin")], by = "curve_id")

##########
###### 2Aiii. Fitting: Sharpe-Schoolfield (6 parameters)  ----
##########

## Since the Sharpe-Schoolfield model has 6 parameters, we need curves with >= 6
## data points from constant temp treatment

curve_needed_sharpe <- num_of_const %>% 
  filter(n >= 6)

curve_needed_sharpe <- c(curve_needed_sharpe$curve_id) 
curve_needed_sharpe # 45 curves with >= 6 const temp treatment

# subset the dataset
data_const_sharpe <- data_const %>% 
  filter(curve_id %in% curve_needed_sharpe)

unique(data_const_sharpe$curve_id) # should be 45 curves


## start progress bar and estimate time it will take
number_of_models <- 1
number_of_curves <- length(unique(data_const_sharpe$curve_id))

## setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models-5, # because of the three curves with errors
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")


## fit the chosen model formulation in rTPC
# Since I got an error with curve_id 11, let's skip it for now
d_fits_sharpe <- nest(data_const_sharpe %>% filter(!is.element(curve_id, c(11,47,42,41,40))), data = c(temp, rate)) %>%
  mutate(sharpeschoolfull = 
           map(data,
               ~nls_multstart_progress(
                 rate ~ sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 20),
                 data = .x,
                 iter = c(3,3,3,3,3,3),  # number of 3's = number of parameters
                 start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
                 start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
                 lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                 upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                 supp_errors = 'Y',
                 convergence_count = FALSE)
               )
         )



##########
###### 2B. Predicting ----
##########

# create new list column of for high resolution data (0.1ºC)
d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(0, 40, by = 0.1)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(briere, quadratic, lactin)) %>% 
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%

  # select only the columns we want to keep
  select(curve_id, study_id, trait, trait_def, units, genus, species, 
         additional_complexity, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(d_preds)


##########
###### 2C. Bootstrapping ----
##########
# get coefs
d_fits <- mutate(d_fits, coefs = map(briere, coef))

# create empty list column
d_fits <- mutate(d_fits, bootstrap = list(rep(NA, n())))

# run for loop to bootstrap each refitted model
for(i in 1:nrow(d_fits)){
  temp_data <- d_fits$data[[i]]
  temp_fit <- nlsLM(rate ~ briere1_1999(temp = temp, tmin, tmax, a),
                    data = temp_data,
                    start = d_fits$coefs[[i]],
                    lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'),
                    upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'))
  boot <- Boot(temp_fit, method = 'residual')
  d_fits$bootstrap[[i]] <- boot
  rm(list = c('temp_fit', 'temp_data', 'boot'))
}

head(d_fits)


# get the raw values of each bootstrap
d_fits <- mutate(d_fits, output_boot = map(bootstrap, function(x) x$t))

# calculate predictions with a gnarly written function
d_fits <- mutate(d_fits, preds = map2(output_boot, data, function(x, y){
  temp <- as.data.frame(x) %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(0, 45, by = 0.5))) %>%
    ungroup() %>%
    mutate(pred = briere1_1999(temp = temp, tmin, tmax, a))
  return(temp)
}))


# select, unnest and calculate 95% CIs of predictions
boot_conf_preds <- select(d_fits, curve_id, preds) %>%
  unnest(preds) %>%
  drop_na(pred) %>% 
  group_by(curve_id, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')


plots_boots <- ggplot() +
  geom_line(data = d_preds, aes(x = temp, y = .fitted, col = model_name)) +
  geom_ribbon(data = boot_conf_preds, aes(temp, ymin = conf_lower, ymax = conf_upper), fill = 'blue', alpha = 0.3) +
  geom_point(data = data_const, aes(temp, rate), size = 2) +
  facet_wrap(~curve_id, scales = 'free_y', ncol = 6) +
  theme_bw(base_size = 12) +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait',
       title = 'All fitted thermal performance curves',
       subtitle = 'Briere in green; quadratic in orange')


# ggsave("figures/tpc-fitted-20250513/tpc-plots-boots.png", plots_boots, width = 45, height = 45, dpi = 300)

##########
###### 2D. Plotting ----
##########

# plot
plots <- ggplot(d_preds) +
  geom_line(aes(x = temp, y = .fitted, col = model_name)) +
  geom_point(data = data_const, aes(temp, rate)) +
  facet_wrap(~curve_id, scales = 'free_y', ncol = 6) +
  theme_bw() +
  #theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait',
       title = 'All fitted thermal performance curves',
       subtitle = 'Briere in green; Lactin in orange, Quadratic in purple')

# ggsave("figures/tpc-fitted-20250513/tpc-plots2.png", plots, width = 45, height = 45, dpi = 300)


## Plot each curve on its own as the above plot is different to read

# Create a directory for output files
output_dir <- "figures/tpc-fitted-20250514"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Generate multiple plots and save them
unique_ids <- unique(d_preds$curve_id)  # Store unique IDs for iteration

for (i in seq_along(unique_ids)) {
  # Create subset of data
  mini_df_pred <- d_preds %>% filter(curve_id == unique_ids[i])
  mini_df_const <- data_const %>% filter(curve_id == unique_ids[i])
  
  # Create a ggplot
  p <- ggplot() +
    geom_line(data = mini_df_pred, aes(x = temp, y = .fitted, col = model_name)) +
    geom_point(data = mini_df_const, aes(temp, rate)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs() +
    xlim(0,40) +
    scale_color_brewer(type = 'qual', palette = 2) +
    labs(title = paste("Curve ID", unique(mini_df_pred$curve_id), "Study ID",
                       unique(mini_df_pred$study_id), unique(mini_df_pred$trait), 
                       unique(mini_df_pred$additional_complexity)),
         x = 'Temperature (ºC)',
         y = paste(unique(mini_df_pred$trait), "(", mini_df_pred$units, ")")) +
    theme_bw()
  
  # Define file path
  file_name <- paste0(output_dir, "/plot_", unique(mini_df_pred$curve_id), "_", 
                      unique(mini_df_pred$study_id), ".png")
  
  # Save the plot
  ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300)
  
  print(paste("Saved:", file_name))
}

