

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(roll)

### make up TPCs

f_curve <- function(x) 0.06*exp(0.09*x)*(1-((x-15)/(34/2))^2)
g_curve <- function(x) 0.01*exp(0.15*x)


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 

tpcs <- p + 
	stat_function(fun = f_curve, color = "purple", size = 2) +
	stat_function(fun = g_curve, color = "orange", size = 2) +
	xlim(0, 35) + ylim(-0.1, 0.5) + 
	xlab("Temperature (°C)") + ylab("Trait") + geom_hline(yintercept = 0)


### make up temperature time series

set.seed(101)
time <- seq_along(1:100)
temps <- data.frame(temperature = runif(100, min = 5, max = 33), time = time)


### plot temps over time

temp_plot <- temps %>% 
	ggplot(aes(x = time, y = temperature)) + geom_line(size = 1) +
	ylab("Temperature (°C)") + xlab("Time (days)")

### get trait values
f_time_series <- data.frame(time = time, ftrait = sapply(X = temps$temperature, FUN = f_curve), temperature = temps$temperature)
g_time_series <- data.frame(time = time, gtrait = sapply(X = temps$temperature, FUN = g_curve), temperature = temps$temperature)


### plot traits over time

both <- left_join(f_time_series, g_time_series)


trait_plot <- both %>% 
	gather(key = trait_type, value = trait_value, ftrait, gtrait) %>% 
	ggplot(aes(x = time, y = trait_value, color = trait_type)) + geom_line(size = 1) +
	scale_color_manual(values = c("purple", "orange")) +ylab("Trait value") + xlab("Time (days)") +
	theme(legend.position = "top")


traits <- both %>% 
	select(ftrait, gtrait) %>% 
	as.matrix()


# rolling covariances
results2 <- roll_cov(traits, width = 2, center = FALSE)
covariances2 <- as.data.frame(results2) %>% 
	select(contains("f")) %>% 
	rownames_to_column(var = "covariance") %>% 
	filter(covariance == "gtrait") %>% 
	gather(key = time_point, value = covariance) %>% 
	mutate(time_point = str_replace_all(time_point, "[a-z, .]", "")) %>% 
	mutate(time_point = as.numeric(time_point)) %>% 
	rename(covariance2 = covariance)

results5 <- roll_cov(traits, width = 5, center = FALSE)
covariances5 <- as.data.frame(results5) %>% 
	select(contains("f")) %>% 
	rownames_to_column(var = "covariance") %>% 
	filter(covariance == "gtrait") %>% 
	gather(key = time_point, value = covariance) %>% 
	mutate(time_point = str_replace_all(time_point, "[a-z, .]", "")) %>% 
	mutate(time_point = as.numeric(time_point)) %>% 
	rename(covariance5 = covariance)

results10 <- roll_cov(traits, width = 10, center = FALSE)
covariances10 <- as.data.frame(results10) %>% 
	select(contains("f")) %>% 
	rownames_to_column(var = "covariance") %>% 
	filter(covariance == "gtrait") %>% 
	gather(key = time_point, value = covariance) %>% 
	mutate(time_point = str_replace_all(time_point, "[a-z, .]", "")) %>% 
	mutate(time_point = as.numeric(time_point)) %>% 
	rename(covariance10 = covariance)

results20 <- roll_cov(traits, width = 20, center = FALSE)
covariances20 <- as.data.frame(results20) %>% 
	select(contains("f")) %>% 
	rownames_to_column(var = "covariance") %>% 
	filter(covariance == "gtrait") %>% 
	gather(key = time_point, value = covariance) %>% 
	mutate(time_point = str_replace_all(time_point, "[a-z, .]", "")) %>% 
	mutate(time_point = as.numeric(time_point)) %>% 
	rename(covariance20 = covariance)



all_cov <- left_join(covariances2, covariances5) %>% 
	left_join(covariances10) %>% 
	left_join(covariances20) %>% 
	gather(key = slice, value = covariance, 2:5) %>% 
	mutate(slice = str_replace_all(slice, "[a-z, .]", "")) %>% 
	rename(time_slice = slice) %>% 
	mutate(time_slice = as.numeric(time_slice))


cov_plot <- all_cov %>% 
	ggplot(aes(x = time_point, y = covariance, color = factor(time_slice), group = time_slice)) + geom_line(size = 1) +
	 scale_color_viridis_d(name = "Duration of time slice") + geom_hline(yintercept = 0) + ylab("cov(ftrait, gtrait)") +
	xlab("Time (days)") + theme(legend.position = "top")


all_plots1 <- plot_grid(temp_plot, trait_plot, cov_plot, nrow = 3, ncol = 1)
all_plots0 <- plot_grid(tpcs, nrow = 1, ncol = 1)

all_plots3 <- plot_grid(all_plots0, all_plots1, rel_widths = c(0.9, 1), rel_heights = c(2, 5), nrow = 1, ncol = 2)
save_plot("figures/covariance-messing.png", all_plots3,  base_height = 7, base_width = 13)
