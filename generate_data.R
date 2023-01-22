generate_AtemEtAl2017 <- function(n, censoring = "light", independent = TRUE, prophaz = TRUE) {
  z <- rbinom(n = n, size = 1, prob = 0.5)
  if (prophaz) {
    x <- rweibull(n = n, shape = 0.75, scale = 0.25)  
  } else {
    x <- rweibull(n = n, shape = 0.6 + 0.9 * z, scale = 0.25)
  }
  e <- rnorm(n = n, mean = 0, sd = 1)
  y <- 1 + 0.5 * x + 0.25 * z + e 
  q <- ifelse(censoring == "light", 2, 0.35)
  c1 <- rweibull(n = n, shape = 0.5, scale = q)
  c2 <- rweibull(n = n, shape = 1, scale = q)
  if (independent) {
    c <- c2
  } else {
    c <- c1
    c[y <= 1.1] <- c2[y <= 1.1]
  }
  w <- pmin(x, c)
  d <- as.numeric(x <= c)
  return(data.frame(z, w, y, d))
}