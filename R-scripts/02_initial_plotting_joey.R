

## Joey's plotting of the TPCs


library(tidyverse)
library(cowplot)


tdata <- read_csv("Data/ExtractedDataAllStudies.csv")


tdata2 <- tdata %>% 
	mutate(study_digits = nchar(study_ID)) %>% 
	filter(study_digits == 3)


tdata2 %>% 
	filter(study_ID == 102) %>% 
	ggplot(aes(x = mean_temp, y = response, color = factor(temp_regime))) + geom_point() +
	facet_wrap( ~ study_ID, scales = "free")
