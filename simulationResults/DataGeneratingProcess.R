source("DGPsetup.R")
source("../estimation/semiparametricAutoregressiveSum.R")


gdp <- function(Kvalue, Tvalue, paramDelta, paramOmega, paramAlpha, paramBeta, vectorSizeIntensity, init, xRates){
   set.seed(1316)
   result <- vector(mode = 'list', length = Kvalue)
   for (kval in 1:Kvalue) {
      resultk <- matrix(data = NA, ncol = 2+length(xRates), nrow = 2*Tvalue)
      resultk[1,1] <- init
      resultk[1, 2] <- 1+rpois(n = 1, lambda = vectorSizeIntensity[[kval]])
      for (tval in 2:(2*Tvalue)) {
         ## Covariates update
         resultk[tval, 3:(2+length(xRates))] <- rexp(n = length(xRates), rate = xRates)
         ## lambda computation
         fregression <- exp(paramOmega + paramAlpha * resultk[tval - 1, 1]/resultk[tval - 1, 2]  +  sum(paramBeta * resultk[tval, 3:(2+length(xRates))]))
         lambda <- log(1 + paramDelta + fregression) 
         ## update n value and y value
         resultk[tval, 2] <-  1+rpois(n = 1, lambda = vectorSizeIntensity[kval])
         resultk[tval, 1] <- sum(rexp(n = resultk[tval, 2], rate = 1/lambda))
      }
      result[[kval]] <- resultk[(Tvalue + 1):(2*Tvalue),]
      result[[kval]] <- as.data.frame(result[[kval]])
      names(result[[kval]]) <- c("yvar", "sizevar", paste0("X", 1:10))
   }
   result
}


## gdp2

gdp2 <- function(Kvalue, Tvalue, paramDelta, paramOmega, paramAlpha, paramBeta, vectorSizeIntensity, init, xRates){
   set.seed(1316)
   result <- vector(mode = 'list', length = Kvalue)
   for (kval in 1:Kvalue) {
      resultk <- matrix(data = NA, ncol = 2+length(xRates), nrow = 2*Tvalue)
      resultk[1,1] <- init
      resultk[1, 2] <- 1+rpois(n = 1, lambda = vectorSizeIntensity[[kval]])
      for (tval in 2:(2*Tvalue)) {
         ## Covariates update
         resultk[tval, 3:(2+length(xRates))] <- rexp(n = length(xRates), rate =  kval * xRates)
         ## lambda computation
         fregression <- exp(paramOmega + paramAlpha * resultk[tval - 1, 1]/resultk[tval - 1, 2]  +  sum(paramBeta * resultk[tval, 3:(2+length(xRates))]))
         lambda <- log(1 + paramDelta + fregression) 
         ## update n value and y value
         resultk[tval, 2] <-  1+rpois(n = 1, lambda = vectorSizeIntensity[kval])
         resultk[tval, 1] <- sum(rexp(n = resultk[tval, 2], rate = 1/lambda))
      }
      result[[kval]] <- resultk[(Tvalue + 1):(2*Tvalue),]
      result[[kval]] <- as.data.frame(result[[kval]])
      names(result[[kval]]) <- c("yvar", "sizevar", paste0("X", 1:10))
   }
   result
}



#### example

sample.1 <- gdp(Kvalue = 1, Tvalue = 200, paramDelta = delta[1], paramAlpha = alpha,
                paramBeta = beta, paramOmega = 1.5, vectorSizeIntensity = intensiteTau[[1]],
                init = initiale, xRates = xCovariatesRates)


## test estimation

x <-  paste0("X", 1:10)

 
model.1 <- modelAutoregressive$new(sample.1)
estimation.1 <- model.1$estimation(y = "yvar", size = "sizevar",
                                   X = x, 
                                   skipto = 2)
erreur.1 <- model.1$standardErrors(mat = FALSE)

 


