## SINGLE IMPUTATION functions
##
## input: missData - a data matrix with NAs for missingData
## output: output - the data matrix with imputed missing values
##X|censored or X|censored, Z
##Kap impute 2 uses density function


kapImpute2 <- function(missData){
  missData <- as.matrix(missData)
  # Count missing values per variable
  ROWna <- rowSums(is.na(missData))
  # Sort row indices by missingness
  sortRowNA <- sort(ROWna, index = TRUE)
  # Find row with minimum missingness
  MissVecStart <- (min(which(sortRowNA$x != 0)))

  nSubj <- dim(missData)[2]
  nVar <- dim(missData)[1]

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  output <- missData
  for (i in MissVecStart:nVar){
    # NOTE: index of ordered minimum count missing is i

    # to get row with respective missing count
    index <- sortRowNA$ix[i]
    rowMissingCount <- sortRowNA$x[i]

    # column indices for NA in row of consideration
    missSubj <- which(is.na(output[index,]))

    # estimating censoring point as min of metabolite
    c <- min(output[index,], na.rm = TRUE)

    for (j in seq_len(length(missSubj))){
      S.est <- 1 - (rowMissingCount/(nVar -rowMissingCount))
      #print(S.est)
      meanEst <- rowMeans(output, na.rm = TRUE)[index]
      sdEst <- sd(output[index,], na.rm = TRUE)

      # E[x_j | x_j < c]
      numerator <- (meanEst/2)*erf((c-meanEst)/(sqrt(2)*sdEst)) + (meanEst/2) -
        (sdEst/sqrt(2*pi))*exp(-((c-meanEst)^2)/(2*sdEst^2))
      denominator <- pnorm(c, meanEst, sdEst)
      #print(denominator)
      output[index, missSubj[j]] <- numerator/denominator

      rowMissingCount = rowMissingCount - 1
    }
  }
  return(output)
}
