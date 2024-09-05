library(ggplot2)
library(wesanderson)

ggplot() + 
  stat_function(fun = dgamma, 
                args = list(shape = 1.5, scale = 2), 
                aes(color = "Gamma"), 
                xlim = c(0, 3)) + 
  stat_function(fun = dlnorm, 
                args = list(meanlog = 0, sdlog = 0.5), 
                aes(color = "Log-Normal"),  
                xlim = c(0, 3)) + 
  stat_function(fun = dnorm, 
                args = list(mean = 1.5, sd = 0.5), 
                aes(color = "Normal"), 
                xlim = c(0, 3))  + 
  stat_function(fun = dweibull, 
                args = list(shape = 0.75, scale = 0.25), 
                aes(color = "Weibull"), 
                xlim = c(0, 3)) + 
  scale_color_manual(name = "Distribution:", 
                     values = c("Gamma" = wes_palette("Zissou1", n = 5, type = "discrete")[1], 
                                "Log-Normal" = wes_palette("Zissou1", n = 5, type = "discrete")[2], 
                                "Normal" = wes_palette("Zissou1", n = 5, type = "discrete")[3], 
                                "t" = wes_palette("Zissou1", n = 5, type = "discrete")[4], 
                                "Weibull" = wes_palette("Zissou1", n = 5, type = "discrete")[5])) + 
  theme_minimal(base_size = 20) + theme(legend.position = "top") + 
  labs(x = "x", y = "Probability Density Function (PDF)")

# Save as 10" wide x 6" tall
ggsave("Figure-DistX.png", width = 5, height = 4, units = "in")
