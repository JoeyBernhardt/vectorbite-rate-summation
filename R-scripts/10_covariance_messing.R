


time_chunk <- function(x){

group <-  sort(rep(seq(from = 1, to = x), times = 100/x))

all <- both %>% 
	mutate(group = group) %>% 
	mutate(slice_length = 100/x) %>% 
	mutate(R_0 = ftrait*gtrait) %>% 
	group_by(group) %>% 
	mutate(`<f(T)g(T)>` = mean(R_0)) %>%
	mutate(`<f(T)>`= mean(ftrait)) %>% 
	mutate(`<g(T)>` = mean(gtrait)) %>% 
	mutate(`Cov[f(T), g(T)]` = cov(ftrait, gtrait)) %>% 
	mutate(`<f(T)><g(T)> + Cov[f(T), g(T)]` = (`<f(T)>` * `<g(T)>`) + `Cov[f(T), g(T)]`) %>% 
	gather(key = estimate, value = value, 8:12)
}


x <- c(1, 2, 4, 5, 10, 20, 25, 50)

all_slices <- x %>% 
	map_df(.f = time_chunk)

all_slices %>% 
	ggplot(aes(x = time, y = value, color = estimate)) + geom_line() +
	facet_wrap( ~ slice_length)

all_slices %>% 
	spread(key = estimate, value = value) %>% 
	mutate(difference = `<f(T)><g(T)> + Cov[f(T), g(T)]` - `<f(T)g(T)>`) %>%
	group_by(slice_length) %>% 
	summarise(average_difference = mean(difference)) %>% 
	ggplot(aes(x  = slice_length, y = average_difference)) + geom_point
