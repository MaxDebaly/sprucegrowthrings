set.seed(3025)
numberOfSitesK <- c(5, 10, 15, 20)
pathLengthT <- c(50, 100, 200)
initiale <- 10
#numberOfCovariates <- 10
omega <- lapply(numberOfSitesK, function(k) seq(from = -0.5*k, to = 0.5*k, length.out = k))
delta <-  0.5
alpha <- 0.6
beta <- c(0, 1, -1, 0.5, -0.5, -1.5, 1.5, -2, 2, 0)
xCovariatesRates <- c(1.2, 2.1, 3.2, 2.3, 4, 4, 2.3, 3.2, 2.1, 1.2)
intensiteTau <-  lapply(numberOfSitesK, function(k) rexp(n = k, rate = 1/k))

