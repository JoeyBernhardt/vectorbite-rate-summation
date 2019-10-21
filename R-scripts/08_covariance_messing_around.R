

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

### testing out the covariance stuff

f_trait <- function(x) 0.05*exp(0.09*x)*(1-((x-15)/(34/2))^2)
g_trait <- function(x) 0.01*exp(0.15*x)


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 

p + 
	stat_function(fun = f_trait, color = "purple") +
	stat_function(fun = g_trait, color = "orange") +
	xlim(0, 35) + ylim(-0.1, 0.5) + 
	theme_bw() + xlab("Temperature (°C)") + ylab("Performance") + geom_hline(yintercept = 0)
