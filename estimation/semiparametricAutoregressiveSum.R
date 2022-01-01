library(R6)

## model for class




likelihoodAutoregressive <- function(listOfData, vectorOfParameters, yvar, ysize,  Xvar, skipto){
   
   nLength <- length(listOfData)
   
   numberOfParameters <- length(vectorOfParameters)
   
   vectorOfLength <- unlist(lapply(listOfData, nrow))
   
   
   individualIntercept <- vectorOfParameters[1:nLength]
   vectorBeta <- vectorOfParameters[(nLength + 1):(nLength+length(Xvar))]
   alphaParameter <- vectorOfParameters[nLength+length(Xvar) + 1]
   deltaParameter <- vectorOfParameters[numberOfParameters]
   
   
   
   result <- 0
   
   for(site in 1:nLength){
      resultSite <- 0
      for(time in skipto:vectorOfLength[site]){
         fregression <- exp(individualIntercept[site]  + alphaParameter * listOfData[[site]][time-1, yvar]/listOfData[[site]][time-1, ysize] + sum(vectorBeta * listOfData[[site]][time,Xvar]))
         lambda <-   log(1+deltaParameter+fregression)
         resultSite <-  resultSite +   listOfData[[site]][time, yvar] / lambda +  listOfData[[site]][time, ysize] * log(lambda)
      }
      result <-  result + resultSite
   }
   
   result
}



stdErrorAutoregressive <- function(listOfData, vectorOfParameters, yvar, ysize, Xvar, skipto){
   
   
   nLength <- length(listOfData)
   
   numberOfParameters <- length(vectorOfParameters)
   
   vectorOfLength <- unlist(lapply(listOfData, nrow))
   
   individualIntercept <- vectorOfParameters[1:nLength]
   vectorBeta <- vectorOfParameters[(nLength + 1):(nLength+length(Xvar))]
   alphaParameter <- vectorOfParameters[nLength+length(Xvar) + 1]
   deltaParameter <- vectorOfParameters[numberOfParameters]
   
   
   I <- 0
   J <- 0
   
   for(site in 1:nLength){
      ISite <- 0
      JSite <- 0
      
      for(time in skipto:vectorOfLength[site]){
         fregression <- exp(individualIntercept[site]  + alphaParameter * listOfData[[site]][time-1, yvar]/listOfData[[site]][time-1, ysize] + sum(vectorBeta * listOfData[[site]][time,Xvar]))
         lambda <-   log(1+deltaParameter+fregression)
         
         partialLambda <- c(unlist(c((1:nLength == site), listOfData[[site]][time,Xvar],  listOfData[[site]][time-1, yvar]/listOfData[[site]][time-1, ysize])) * (fregression /(1+ deltaParameter + fregression)), 1/(1+ deltaParameter + fregression))
         smallS <- (1/lambda) * partialLambda  * (
            listOfData[[site]][time,ysize] -(listOfData[[site]][time,yvar] / lambda)
         )
         smallH <- (listOfData[[site]][time,ysize]/(lambda^2)) * partialLambda %*% t(partialLambda) 
         
         ISite <- ISite + smallS %*% t(smallS)
         JSite <- JSite + smallH
      }
      
      I <- I + ISite/(vectorOfLength[site]-skipto + 1)
      J <- J + JSite/(vectorOfLength[site]-skipto + 1)
   }
   
   return(list(Imatrix = I, Jmatrix = J))
}



modelAutoregressive <- R6Class(classname = "modelAutoregressive",
                                            public = list(
                                               
                                               dataList = NA,
                                               parameters = NA,
                                               stderrors = NA,
                                               numberOfParameters = NA,
                                               yvar = NA,
                                               Xvar = NA,
                                               ysize = NA,
                                               skipto = NA,
                                               
                                               
                                               initialize = function(data){
                                                  self$dataList <-  data
                                                  
                                               },
                                               
                                               estimation = function(initialParameters = NULL, y, size,  X, skipto){
                                                  self$yvar <-  y
                                                  self$Xvar <- X
                                                  self$ysize <- size
                                                  self$skipto <- skipto
                                                  lenBlock <-  length(self$dataList)
                                                  if(is.null(initialParameters)) initialParameters = c(rep(0, lenBlock), rep(0, length(X)),  runif(2))
                                                  
                                                  optimization <- optim(par = initialParameters, fn = likelihoodAutoregressive, listOfData = self$dataList, 
                                                                        yvar = y, ysize = size, skipto = skipto,  Xvar = X, lower = c(rep(-Inf, lenBlock+length(X) + 1), 0),
                                                                        method = "L-BFGS-B")
                                                  
                                                  self$parameters <- optimization$par
                                                  return(self$parameters)
                                               },
                                               
                                               standardErrors = function(mat = FALSE){
                                                  matrices = stdErrorAutoregressive(listOfData = self$dataList, 
                                                                                    vectorOfParameters = self$parameters, 
                                                                                    yvar = self$yvar, ysize = self$ysize,
                                                                                    Xvar = self$Xvar, skipto = self$skipto)
                                                  Imatrix <-  matrices$Imatrix
                                                  invJmatrix <- solve(matrices$Jmatrix)
                                                  self$stderrors <-  sqrt(diag(invJmatrix %*% Imatrix %*% t(invJmatrix)/(nrow(self$dataList[[1]]) - self$skipto + 1)))
                                                  
                                                  if(mat){
                                                     return(list(matrices, self$stderrors))
                                                  }
                                                  else{
                                                     return(self$stderrors)
                                                  }
                                               }
                                            )
) 



