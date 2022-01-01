source("DGPsetup.R")
source("../estimation/semiparametricAutoregressiveSum.R")

indexesSimulation <- expand.grid(K = numberOfSitesK, Time = pathLengthT)
indexesSimulation <- indexesSimulation[1:11,]

## simulation 1
resultsSimulation <- apply(indexesSimulation, 1, function(indexes) {
   print(indexes)
   sample <- gdp(indexes[1],indexes[2],
                 paramDelta = delta[1], paramAlpha = alpha,
                 paramBeta = beta, paramOmega = omega[which(numberOfSitesK == indexes[1])][[1]], vectorSizeIntensity = intensiteTau[[which(numberOfSitesK == indexes[1])]],
                 init = initiale, xRates =  xCovariatesRates)
   model <- modelAutoregressive$new(sample)
   estimation <- model$estimation(y = "yvar", size = "sizevar",
                                  X = x, 
                                  skipto = 2)
   erreur <- model$standardErrors(mat = FALSE)
   return(rbind(estimation, erreur))
})


#result 1
save(resultsSimulation, file = "resultsSimulation.RData")


## simulation 2
resultsSimulation2 <- apply(indexesSimulation, 1, function(indexes) {
   print(indexes)
   sample <- gdp2(indexes[1],indexes[2],
                  paramDelta = delta[1], paramAlpha = alpha,
                  paramBeta = beta, paramOmega = omega[which(numberOfSitesK == indexes[1])][[1]], vectorSizeIntensity = intensiteTau[[which(numberOfSitesK == indexes[1])]],
                  init = initiale, xRates =  xCovariatesRates)
   model <- modelAutoregressive$new(sample)
   estimation <- model$estimation(y = "yvar", size = "sizevar",
                                  X = x, 
                                  skipto = 2)
   erreur <- model$standardErrors(mat = FALSE)
   return(rbind(estimation, erreur))
})


 

## resultat 2

save(resultsSimulation2, file = "resultsSimulation2.RData")

