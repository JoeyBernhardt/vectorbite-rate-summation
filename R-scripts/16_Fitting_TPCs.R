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
##    3. Model Selection & Model Averaging
##    4. Bootstrapping
##    5. Plotting



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
library(ggrepel)

setwd("~/Documents/UofG/vectorbite-rate-summation")

# load in data
data <- read_csv("data-processed/vb_maggie_data_20250611.csv")

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

## Since the Briere and Quadratic model has 3 parameters, 
## we need curves with >= 4 data points from constant temp treatment
num_of_const <- data %>% 
  filter(dtr == 0) %>% 
  group_by(curve_id) %>% 
  summarise(n = n()) %>% 
  arrange(n)

curve_needed <- num_of_const %>% 
  filter(n >= 4)

curve_needed <- c(curve_needed$curve_id) 
curve_needed # 60 curves with >= 4 const temp treatment


data_const <- data_const %>% 
  filter(curve_id %in% curve_needed)

unique(data_const$curve_id) # should be 60 curves


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
# d_fits_bri_quad <- nest(data_const, data = c(temp, rate)) %>%
#   mutate(briere = map(data, 
#                       ~nls_multstart_progress(
#                         rate~briere1_1999(temp = temp, tmin, tmax, a),
#                         data = .x,
#                         iter = c(3,3,3),
#                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere1_1999') - 10,
#                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere1_1999') + 10,
#                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere1_1999'),
#                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere1_1999'),
#                         supp_errors = 'Y',
#                         convergence_count = FALSE)
#                       ), 
#          quadratic = map(data, 
#                       ~nls_multstart_progress(
#                         formula = rlang::expr(rate ~ quadratic(temp, tmin, tmax, a)), #rlang::expr is very very important. Without it the code below won't run
#                         data = .x,
#                         iter = c(3,3,3),
#                         start_lower = quadratic.starting_vals(.x$temp, .x$rate) - 10,
#                         start_upper = quadratic.starting_vals(.x$temp, .x$rate) + 10,
#                         lower = quadratic.lower_lims(.x$temp, .x$rate),
#                         upper = quadratic.upper_lims(.x$temp, .x$rate),
#                         supp_errors = 'Y',
#                         convergence_count = FALSE)
#                       )
#          )


# save(d_fits_bri_quad, file = "R-scripts/model-output/d_fits_bri_quad.Rdata")



##########
###### 2Aii. Fitting: Lactin (4 parameters)  ----
##########

## Since the Lactin model has 4 parameters, we need curves with >= 5 data points from constant temp treatment
curve_needed_lac <- num_of_const %>% 
  filter(n >= 5)

curve_needed_lac <- c(curve_needed_lac$curve_id) 
curve_needed_lac # 54 curves with >= 5 const temp treatment


data_const_lac <- data_const %>% 
  filter(curve_id %in% curve_needed_lac)

unique(data_const_lac$curve_id) # should be 54 curves


## start progress bar and estimate time it will take
number_of_models <- 1
number_of_curves <- length(unique(data_const_lac$curve_id))

## setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")


## fit the chosen model formulation in rTPC
# d_fits_lac <- nest(data_const_lac, data = c(temp, rate)) %>%
#   mutate(lactin = map(data, 
#                       ~nls_multstart_progress(
#                         rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
#                         data = .x,
#                         iter = c(3,3,3,3),
#                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
#                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
#                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
#                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
#                         supp_errors = 'Y',
#                         convergence_count = FALSE)
#                       )
#          )

# save(d_fits_lac, file = "R-scripts/model-output/d_fits_lac.Rdata")


##########
###### 2Aiii. Fitting: Sharpe-Schoolfield (6 parameters)  ----
##########

## Since the Sharpe-Schoolfield model has 6 parameters, we need curves with >= 7

## Data points from constant temp treatment
curve_needed_sharpe <- num_of_const %>% 
  filter(n >= 7)

curve_needed_sharpe <- c(curve_needed_sharpe$curve_id) 
curve_needed_sharpe # 36 curves with >= 7 const temp treatment

## subset the dataset
data_const_sharpe <- data_const %>% 
  filter(curve_id %in% curve_needed_sharpe)

unique(data_const_sharpe$curve_id) # should be 36 curves


## Get all the curve_ids
curves_id_sharpe <- unique(data_const_sharpe$curve_id)


## Create nested dataframe like the format we used above
d_fits_sharpe <- nest(data_const_sharpe, data = c(temp, rate))

## Initialize a list-column to hold  model fits 
d_fits_sharpe$sharpeschoolfull <- vector("list", nrow(d_fits_sharpe))


## We will fit Sharpe-Schoolfield model one-by-one using for loops (since the map function hates Sharpes-Schoolfield)
## curve_id 9, 48, and 67 will be skipped in the following code since they do not have eh starting value;

# for (id in curves_id_sharpe) {
# 	print(paste("Fitting curve_id", id))
# 	
# 	# Extract the data of each curve to fit the model
# 	data <- subset(data_const_sharpe, curve_id == id)
# 	
# 	# get start values and fit model
# 	start_vals <- get_start_vals(data$temp, data$rate, model_name = 'sharpeschoolfull_1981')
# 	
# 	## Check for NA in starting values
# 	if (any(is.na(start_vals))) {
# 		message("  -> skipping curve_id ", id, " due to NA in starting values")
# 		next # skip to next curve
# 	}
# 	
# 	# fit model, use tryCatch to handle failures
# 	fit <- tryCatch(
# 		{
# 			nls_multstart(
# 				rate~sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 20),
# 				data = data,
# 				iter = c(3,3,3,3,3,3),
# 				start_lower = start_vals - 10,
# 				start_upper = start_vals + 10,
# 				lower = get_lower_lims(data$temp, data$rate, model_name = 'sharpeschoolfull_1981'),
# 				upper = get_upper_lims(data$temp, data$rate, model_name = 'sharpeschoolfull_1981'),
# 				supp_errors = 'Y',
# 				convergence_count = FALSE
# 				)
# 		},
# 		error = function(e) {
# 			message("  -> fitting failed for curve_id ", id, ": ", e$message)
# 			return(NA)
# 		}
# 	)
# 	
# 	# Store the fit if successful
# 	d_fits_sharpe$sharpeschoolfull[d_fits_sharpe$curve_id == id] <- list(fit)
# 	
# }


# save(d_fits_sharpe, file = "R-scripts/model-output/d_fits_sharpe.Rdata")


## Load all the model fits
load("R-scripts/model-output/d_fits_bri_quad.Rdata")
load("R-scripts/model-output/d_fits_lac.Rdata")
load("R-scripts/model-output/d_fits_sharpe.Rdata")


## Combine model fits from the Briere, Quadratic, Lactin and Sharpe-Schoolfield models
d_fits <- d_fits_bri_quad %>% 
	full_join(d_fits_lac[,c("curve_id", "lactin")],  by = "curve_id") %>% 
	full_join(d_fits_sharpe[,c("curve_id", "sharpeschoolfull")],  by = "curve_id")



##########
###### 2B. Predicting ----
##########
temp_grad <- round(seq(0, 50, 0.1), 1) 

# create new list column of for high resolution data (0.1ºC)
d_stack <- mutate(d_fits, new_data = map(data, ~tibble(temp = temp_grad))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(briere, quadratic, lactin, sharpeschoolfull)) 


d_preds <- d_stack %>% 
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%

  # select only the columns we want to keep
  select(curve_id, study_id, trait, trait_def, units, genus, species, 
         additional_complexity, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)



# save(d_preds, file = "R-scripts/model-output/d_preds.Rdata")

##########
###### 3. Model Selection & Model Averaging ----
##########

## Get measures of relative model fit
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance), #AIC and BIC
         AICc = map_dbl(fit, ~ ifelse(!is.null(.x), MuMIn::AICc(.x), NA))) %>% #AICc
  select(-fit) %>%
  unnest(info, keep_empty = T) %>%
  select(curve_id, study_id, trait, model_name, sigma, isConv, AIC, AICc, BIC, df.residual)


            
## We use AICc score to compare between models
## The model with the lowest AICc score is chosen as the model that best supports the data
best_model <- d_ic %>% 
  filter(!is.na(AICc)) %>% 
  filter(isConv == T) %>% # to remove models that do not converge
  group_by(curve_id) %>% 
  slice_min(AICc, with_ties = F) %>% # select the models with the lowest AICc
  ungroup() %>% 
  select(curve_id, study_id, trait, model_name, AICc)

best_model


## Select predicted values of the best model for each curve
d_preds_best <- d_preds %>% 
  semi_join(best_model, by = c("curve_id", "model_name"))

# save(d_preds_best, file = "R-scripts/model-output/d_preds_best.Rdata")


##########
###### 4. Bootstrapping ----
##########

## get coefs
d_fits_bri <- mutate(d_fits_bri_quad[, -11] , coefs = map(briere, coef))
d_fits_quad <- mutate(d_fits_bri_quad[, -10], coefs = map(quadratic, coef))
d_fits_lac <- mutate(d_fits_lac, coefs = map(lactin, coef))
d_fits_sharpe <- mutate(d_fits_sharpe, coefs = map(sharpeschoolfull, coef))


## create empty list column
d_fits_bri <- mutate(d_fits_bri, bootstrap = list(rep(NA, n())))
d_fits_quad <- mutate(d_fits_quad, bootstrap = list(rep(NA, n())))
d_fits_lac <- mutate(d_fits_lac, bootstrap = list(rep(NA, n())))
d_fits_sharpe <- mutate(d_fits_sharpe, bootstrap = list(rep(NA, n())))


## run for loop to bootstrapping Briere model ----
for(i in 1:nrow(d_fits_bri)){
  temp_data <- d_fits_bri$data[[i]]
  temp_fit <- nlsLM(rate ~ briere1_1999(temp = temp, tmin, tmax, a),
                    data = temp_data,
                    start = d_fits_bri$coefs[[i]],
                    lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'),
                    upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'briere1_1999'))
  boot <- Boot(temp_fit, method = 'residual')
  d_fits_bri$bootstrap[[i]] <- boot
  rm(list = c('temp_fit', 'temp_data', 'boot'))
}

head(d_fits_bri)


## run for loop to bootstrapping quadratic model ----
for(i in 1:nrow(d_fits_quad)){
	temp_data <- d_fits_quad$data[[i]]
	temp_fit <- nlsLM(formula = rlang::expr(rate ~ quadratic(temp, tmin, tmax, a)),
					  data = temp_data,
					  start = d_fits_quad$coefs[[i]],
					  lower = quadratic.lower_lims(temp_data$temp, temp_data$rate),
					  upper = quadratic.upper_lims(temp_data$temp, temp_data$rate))
	boot <- Boot(temp_fit, method = 'residual')
	d_fits_quad$bootstrap[[i]] <- boot
	rm(list = c('temp_fit', 'temp_data', 'boot'))
}

head(d_fits_quad)


## run for loop to bootstrapping Lactin model ----
for(i in 1:nrow(d_fits_lac)){
	temp_data <- d_fits_lac$data[[i]]
	temp_fit <- nlsLM(rate ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
					  data = temp_data,
					  start = d_fits_lac$coefs[[i]],
					  lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'lactin2_1995'),
					  upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'lactin2_1995'))
	boot <- Boot(temp_fit, method = 'residual')
	d_fits_lac$bootstrap[[i]] <- boot
	rm(list = c('temp_fit', 'temp_data', 'boot'))
}

head(d_fits_lac)



## run for loop to bootstrapping Sharpe-Schoolfield model ----
## First remove rows that did not fit the Sharpe-Schoolfield model due to missing starting value
d_fits_sharpe <- d_fits_sharpe %>% 
	filter(map_lgl(sharpeschoolfull, ~!is.null(.x)))

for(i in 1:nrow(d_fits_sharpe)){
	temp_data <- d_fits_sharpe$data[[i]]
	temp_fit <- nlsLM(rate ~ sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 20),
					  data = temp_data,
					  start = d_fits_sharpe$coefs[[i]],
					  lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'sharpeschoolfull_1981'),
					  upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'sharpeschoolfull_1981'))
	boot <- Boot(temp_fit, method = 'residual')
	d_fits_sharpe$bootstrap[[i]] <- boot
	rm(list = c('temp_fit', 'temp_data', 'boot'))
}

head(d_fits_sharpe)




## Get the raw values of each bootstrap
d_fits_bri <- mutate(d_fits_bri, output_boot = map(bootstrap, function(x) x$t))
d_fits_quad <- mutate(d_fits_quad, output_boot = map(bootstrap, function(x) x$t))
d_fits_lac <- mutate(d_fits_lac, output_boot = map(bootstrap, function(x) x$t))
d_fits_sharpe <- mutate(d_fits_sharpe, output_boot = map(bootstrap, function(x) x$t))


## Calculate predictions with a gnarly written function
d_fits_bri <- mutate(d_fits_bri, preds = map2(output_boot, data, function(x, y){
  temp <- as.data.frame(x) %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = temp_grad)) %>%
    ungroup() %>%
    mutate(pred = briere1_1999(temp = temp, tmin, tmax, a))
  return(temp)
}))


d_fits_quad <- mutate(d_fits_quad, preds = map2(output_boot, data, function(x, y){
	temp <- as.data.frame(x) %>%
		drop_na() %>%
		mutate(iter = 1:n()) %>%
		group_by_all() %>%
		do(data.frame(temp = temp_grad)) %>%
		ungroup() %>%
		mutate(pred = quadratic(temp, tmin, tmax, a))
	return(temp)
}))


d_fits_lac <- mutate(d_fits_lac, preds = map2(output_boot, data, function(x, y){
	temp <- as.data.frame(x) %>%
		drop_na() %>%
		mutate(iter = 1:n()) %>%
		group_by_all() %>%
		do(data.frame(temp = temp_grad)) %>%
		ungroup() %>%
		mutate(pred = lactin2_1995(temp = temp, a, b, tmax, delta_t))
	return(temp)
}))

d_fits_sharpe <- mutate(d_fits_sharpe, preds = map2(output_boot, data, function(x, y){
	temp <- as.data.frame(x) %>%
		drop_na() %>%
		mutate(iter = 1:n()) %>%
		group_by_all() %>%
		do(data.frame(temp = temp_grad)) %>%
		ungroup() %>%
		mutate(pred = sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 20))
	return(temp)
}))


# select, unnest and calculate 95% CIs of predictions
boot_conf_preds_bri <- select(d_fits_bri, curve_id, preds) %>%
  unnest(preds) %>%
  drop_na(pred) %>% 
  group_by(curve_id, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop') %>% 
	mutate(model_name = "briere")


boot_conf_preds_quad <- select(d_fits_quad, curve_id, preds) %>%
	unnest(preds) %>%
	drop_na(pred) %>% 
	group_by(curve_id, temp) %>%
	summarise(conf_lower = quantile(pred, 0.025),
			  conf_upper = quantile(pred, 0.975),
			  .groups = 'drop') %>% 
	mutate(model_name = "quadratic")


boot_conf_preds_lac <- select(d_fits_lac, curve_id, preds) %>%
	unnest(preds) %>%
	drop_na(pred) %>% 
	group_by(curve_id, temp) %>%
	summarise(conf_lower = quantile(pred, 0.025),
			  conf_upper = quantile(pred, 0.975),
			  .groups = 'drop') %>% 
	mutate(model_name = "lactin")


boot_conf_preds_sharpe <- select(d_fits_sharpe, curve_id, preds) %>%
	unnest(preds) %>%
	drop_na(pred) %>% 
	group_by(curve_id, temp) %>%
	summarise(conf_lower = quantile(pred, 0.025),
			  conf_upper = quantile(pred, 0.975),
			  .groups = 'drop')%>% 
	mutate(model_name = "sharpeschoolfull")


## Now put everything into a dataframe
boot_conf_preds <- rbind(boot_conf_preds_bri, boot_conf_preds_quad, boot_conf_preds_lac, boot_conf_preds_sharpe)


# save(d_fits_bri, file = "R-scripts/model-output/d_fits_bri.Rdata")
# save(d_fits_quad, file = "R-scripts/model-output/d_fits_quad.Rdata")
# save(d_fits_lac, file = "R-scripts/model-output/d_fits_lac.Rdata")
# save(d_fits_sharpe, file = "R-scripts/model-output/d_fits_sharpe.Rdata")
# save(boot_conf_preds, file = "R-scripts/model-output/boot_conf_preds.Rdata")

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


# ggsave("figures/Lilian/tpc-fitted-20250618/tpc-plots-boots.png", plots_boots, width = 45, height = 45, dpi = 300)


## Select predicted values and confidence interval of the best model for each curve
boot_conf_preds_best <- boot_conf_preds %>% 
	semi_join(best_model, by = c("curve_id", "model_name"))

d_preds_best <- d_preds_best %>% 
	inner_join(boot_conf_preds_best, by = c("curve_id", "temp", "model_name"))

# save(d_preds_best, file = "R-scripts/model-output/d_preds_best.Rdata")

##########
###### 5. Plotting ----
##########

# Define custom colors for each model
model_colors <- c(
	"briere" = "#999999",
	"quadratic" = "#E69F00",
	"lactin" = "#56B4E9",
	"sharpeschoolfull" = "#009E73")

# Define custom labels
model_labels <- c(
	"briere" = "Briere",
	"quadratic" = "Quadratic",
	"lactin" = "Lactin II",
	"sharpeschoolfull" = "Sharpe-Schoolfield (full)")


# plot all 4 TPCs of each curve in a big plot ----
plots <- ggplot(d_preds) +
	geom_line(aes(x = temp, y = .fitted, col = model_name)) +
	geom_point(data = data_const, aes(temp, rate)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	facet_wrap(~curve_id, scales = 'free_y', ncol = 6) +
	theme_bw() +
	scale_color_manual(name = "Model", values = model_colors, labels = model_labels) +
	scale_fill_manual(name = "Model", values = model_colors, labels = model_labels) +	
    labs(x = 'Temperature (ºC)',
    	 y = 'Trait',
         title = 'All fitted thermal performance curves')

# ggsave("figures/Lilian/tpc-fitted-20250618/tpc-plots.png", plots, width = 45, height = 45, dpi = 300)


## Plot each curve on its own as the above plot is difficult to read
# Create a directory for output files
output_dir <- "figures/Lilian/tpc-fitted-20250618"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


## Generate multiple plots and save them (all models and without confidence interval)
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
    ylim(0, NA) +
  	scale_color_manual(name = "Model", values = model_colors, labels = model_labels) +
  	scale_fill_manual(name = "Model", values = model_colors, labels = model_labels) +
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


## plot the best fitting TPC of each curve ----
plots <- ggplot(d_preds_best) +
	geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = model_name), alpha = 0.3) +
	geom_line(aes(x = temp, y = .fitted, col = model_name)) +
	geom_point(data = data_const, aes(temp, rate)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	facet_wrap(~curve_id, scales = 'free_y', ncol = 6) +
	theme_bw() +
	scale_color_manual(name = "Model", values = model_colors, labels = model_labels) +
	scale_fill_manual(name = "Model", values = model_colors, labels = model_labels) +	
	labs(x = 'Temperature (ºC)',
		 y = 'Trait',
		 title = 'All fitted thermal performance curves')

# ggsave("figures/Lilian/tpc-fitted-20250618/tpc-plots-best-mod.png", plots, width = 45, height = 45, dpi = 300)



## Plot each curve on its own as the above plot is difficult to read

# Create a directory for output files
output_dir <- "figures/Lilian/tpc-fitted-20250618"

if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}


unique_ids <- unique(d_preds$curve_id)  # Store unique IDs for iteration

for (i in seq_along(unique_ids)) {
	# Create subset of data
	mini_df_const <- data_const %>% filter(curve_id == unique_ids[i])
	mini_df_pred <- d_preds_best %>% filter(curve_id == unique_ids[i])
	mini_df_conf <- boot_conf_preds_best %>% filter(curve_id == unique_ids[i])
	
	# Create a ggplot
	p <- ggplot() +
		geom_line(data = mini_df_pred, aes(x = temp, y = .fitted, col = model_name)) +
		geom_ribbon(data = mini_df_conf, aes(temp, ymin = conf_lower, ymax = conf_upper, fill = model_name), alpha = 0.3) +
		geom_point(data = mini_df_const, aes(temp, rate)) +
		geom_hline(yintercept = 0, linetype = "dashed") +
		ylim(0, NA) +
		scale_color_manual(name = "Model", values = model_colors, labels = model_labels) +
		scale_fill_manual(name = "Model", values = model_colors, labels = model_labels) +
		labs(title = paste("Curve ID", unique(mini_df_pred$curve_id), "Study ID",
						   unique(mini_df_pred$study_id), unique(mini_df_pred$trait), 
						   unique(mini_df_pred$additional_complexity)),
			 x = 'Temperature (ºC)',
			 y = paste(unique(mini_df_pred$trait), "(", mini_df_pred$units, ")")) +
		theme_bw()
	
	# Define file path
	file_name <- paste0(output_dir, "/plot_", unique(mini_df_pred$curve_id), "_", 
						unique(mini_df_pred$study_id), "_best.png")
	
	# Save the plot
	ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300)
	
	print(paste("Saved:", file_name))
}

## Plot curve_id 33 and 34 again (since the Lactin model prediction is very very high at high temp)

