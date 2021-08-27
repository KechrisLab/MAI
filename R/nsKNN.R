nsKNN = function(data_miss,
                 k,
                 iters=1,
                 weighted = FALSE,
                 scale = FALSE,
                 shuffle = TRUE){
  if (scale){
    data_miss_scale = scale(data_miss)
    data_miss = data_miss_scale
  }


  findKneighbors = function(distances, k, col){
    distances = distances[distances[,1]==col | distances[,3]==col, ]
    sortedDistances = distances[order(distances[,2]),]
    neighbors = sortedDistances[seq_len(k), ]
    n = list()
    for(i in seq_len(nrow(neighbors))){n[i] = ifelse(neighbors[i,1] == col,
                                              list(neighbors[i, c(3,2)]),
                                              list(neighbors[i, c(1,2)]))}
    n = matrix(unlist(n), nrow = k, ncol = 2, byrow = TRUE)
    return(n)
  }



  distances = as.matrix(dist(t(data_miss))/sqrt(nrow(data_miss)))
  distances = distances[lower.tri(distances)]
  distances = data.frame(col1 = combn(paste(seq_len(ncol(data_miss)),
                                            sep = ""), 2)[1,],
                         distance = distances,
                         col2 = combn(paste(seq_len(ncol(data_miss)),
                                            sep = ""), 2)[2,])
  distances = apply(distances, 2, as.numeric)

  minVals = apply(data_miss, 1, min, na.rm=TRUE)

  data_imputed = data_miss

  for (col in seq_len(ncol(data_miss))){
    neighbors = findKneighbors(distances, k, col)
    missingNum = which(is.na(data_miss[,col]))

    if(length(missingNum)==0){
      next
    }

    neighborVals = data_miss[missingNum, c(col, neighbors[,1])]
    mins = minVals[missingNum]
    if (is.null(dim(neighborVals))){
      listns = list(c(neighborVals[-1], mins))
      listns = lapply(listns, function(x){ifelse(is.na(x), x[k+1], x)})

      filledNeighborVals = unlist(listns)[-length(neighborVals)]
    } else {
      listns = asplit(cbind(neighborVals[,-1], mins), 1)
      listns = lapply(listns, function(x){ifelse(is.na(x), x[k+1], x)})


    filledNeighborVals = cbind(matrix(unlist(listns),
                                      ncol = ncol(neighborVals),
                                      nrow = nrow(neighborVals),
                                      byrow = TRUE))[,-ncol(neighborVals)]
    }

    weightMultiplier = (1/neighbors[,2])/sum(1/neighbors[,2])

    if (weighted){
      if (is.null(dim(neighborVals))){
        data_imputed[missingNum, col] = sum(
          filledNeighborVals%*%diag(weightMultiplier),
          na.rm=TRUE)
      } else {
        data_imputed[missingNum, col] = apply(
          filledNeighborVals%*%diag(weightMultiplier),
          1, sum, na.rm=TRUE)
      }

    } else {
      if (is.null(dim(neighborVals))){
        data_imputed[missingNum, col] = mean(
          filledNeighborVals%*%diag(weightMultiplier),
          na.rm=TRUE)
      }else{
        data_imputed[missingNum, col] = apply(filledNeighborVals,
                                              1, mean, na.rm=TRUE)
      }

    }

  }

  if (iters < 2){
    if (scale){
      data_imputed = t(apply(data_imputed, 1, function(r){
        r*attr(data_miss_scale,'scaled:scale') +
          attr(data_miss_scale, 'scaled:center')
        }))
      return(data_imputed)
    } else {
      return(data_imputed)
    }

  } else {
    for (j in seq_len((iters-1))){
      if (shuffle){
        numCols = sample(seq_len(ncol(data_imputed)),
                         ncol(data_imputed), replace = FALSE)
      } else {
        numCols = seq_len(ncol(data_imputed))
      }

      for (numCol in numCols){
        neighborCols = numCols[numCols!=numCol]
        distances = as.matrix(dist(t(data_imputed[,c(numCol, neighborCols)])))
        distances = cbind(neighborCols, distances[1,][-1])
        sortDistances = distances[order(distances[,2]),]
        NN = sortDistances[seq_len(k)]
        missingNum = which(is.na(data_miss[,numCol]))
        neighborVals = data_imputed[missingNum, c(NN)]
        weightMultiplier = (1/sortDistances[,2][seq_len(k)])/sum(1/sortDistances[,2][seq_len(k)])
        if (weighted){
          if (is.null(dim(neighborVals))){
            data_imputed[missingNum, numCol] = sum(
              as.matrix(neighborVals)%*%diag(weightMultiplier),
              na.rm=TRUE)
          } else {
            data_imputed[missingNum, numCol] = apply(
              as.matrix(neighborVals)%*%diag(weightMultiplier),
              1, sum, na.rm=TRUE)
          }

        } else {
          if (is.null(dim(neighborVals))){
            data_imputed[missingNum, numCol] = mean(neighborVals,
                                                    na.rm=TRUE)
          } else {
            data_imputed[missingNum, numCol] = apply(neighborVals,
                                                     1, mean, na.rm=TRUE)
          }

        }

      }

    }
   if (scale){
     data_imputed = t(apply(data_imputed,
                            1, function(r){
                              r*attr(data_miss_scale,'scaled:scale') +
                                attr(data_miss_scale, 'scaled:center')
                              }))
     return(data_imputed)
   } else {
     return(data_imputed)
   }

  }

}
