hist_with_density = function(data, func, start = NULL){
  # load libraries
  library(VGAM); library(fitdistrplus); library(ggplot2)
  
  # fit density to data
  fit   = fitdist(data, func, start = start)
  args  = as.list(fit$estimate)
  dfunc = match.fun(paste('d', func, sep = ''))
  
  # plot histogram, empirical and fitted densitiesdat
  p0 = qplot(data, geom = 'blank') +
    geom_line(aes(y = ..density..,colour = 'Empirical'),stat = 'density') +
    stat_function(fun = dfunc, args = args, aes(colour = func))  +
    geom_histogram(aes(y = ..density..), alpha = 0.4) +
    scale_colour_manual(name = '', values = c('red', 'blue')) + 
    theme(legend.position = 'top', legend.direction = 'horizontal')
  return(p0)  
}