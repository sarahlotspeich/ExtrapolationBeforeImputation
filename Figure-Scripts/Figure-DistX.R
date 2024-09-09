library(ggplot2)
library(wesanderson)

# Plot PDFs
plot_pdf = ggplot() + 
  stat_function(fun = dgamma, 
                args = list(shape = 1.5, scale = 2), 
                aes(color = "Gamma"), 
                xlim = c(0, 5)) + 
  stat_function(fun = dlnorm, 
                args = list(meanlog = 0, sdlog = 0.5), 
                aes(color = "Log-Normal"),  
                xlim = c(0, 5)) + 
  stat_function(fun = dweibull, 
                args = list(shape = 0.75, scale = 0.25), 
                aes(color = "Weibull"), 
                xlim = c(0, 5)) + 
  scale_color_manual(name = "Distribution for X:", 
                     values = c("Gamma" = wes_palette("Zissou1", n = 3, type = "discrete")[1], 
                                "Log-Normal" = wes_palette("Zissou1", n = 3, type = "discrete")[3], 
                                "Weibull" = wes_palette("Zissou1", n = 5, type = "discrete")[5])) + 
  theme_minimal(base_size = 20) + theme(legend.position = "top") + 
  labs(x = "x", y = "Probability Density Function")

# Plot hazard functions
hgamma = function(x, shape, scale) {
  dgamma(x = x, shape = shape, scale = scale) / 
    pgamma(q = x, shape = shape, scale = scale, lower.tail = FALSE)
}

hlnorm = function(x, meanlog, sdlog) {
  dlnorm(x = x, meanlog = meanlog, sdlog = sdlog) / 
    plnorm(q = x, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE)
}

hweibull = function(x, shape, scale) {
  dweibull(x = x, shape = shape, scale = scale) / 
    pweibull(q = x, shape = shape, scale = scale, lower.tail = FALSE)
}

plot_haz = ggplot() + 
  stat_function(fun = hgamma, 
                args = list(shape = 1.5, scale = 2), 
                aes(color = "Gamma"), 
                xlim = c(0, 5)) + 
  stat_function(fun = hlnorm, 
                args = list(meanlog = 0, sdlog = 0.5), 
                aes(color = "Log-Normal"),  
                xlim = c(0, 5)) + 
  stat_function(fun = hweibull, 
                args = list(shape = 0.75, scale = 0.25), 
                aes(color = "Weibull"), 
                xlim = c(0, 5)) + 
  scale_color_manual(name = "Distribution for X:", 
                     values = c("Gamma" = wes_palette("Zissou1", n = 3, type = "discrete")[1], 
                                "Log-Normal" = wes_palette("Zissou1", n = 3, type = "discrete")[3], 
                                "Weibull" = wes_palette("Zissou1", n = 5, type = "discrete")[5])) + 
  theme_minimal(base_size = 20) + theme(legend.position = "top") + 
  labs(x = "x", y = "Hazard Function") 

# Combine them and save
ggpubr::ggarrange(plot_pdf, plot_haz, ncol = 1, labels = "AUTO", common.legend = TRUE)

# Save as 10" wide x 6" tall
ggsave("Figure-DistX.png", width = 10, height = 10, units = "in")
