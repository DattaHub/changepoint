n <- 497
theta <- numeric(n)
theta[1:138] <- -0.18
theta[139:225] <- 0.08
theta[226:242] <- 1.07
theta[243:299] <- -0.53
theta[300:308] <- 0.16
theta[309:332] <- -0.69
theta[333:n] <- -0.16
n <- length(theta)
sig2 = 0.04
y <- theta + sqrt(sig2) * rnorm(n)
plot(y, cex=0.5, col="gray")
lines(theta)
