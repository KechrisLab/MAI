# Function to estimate thresholds
check_distance = function(data_miss, data_sub, PercentMiss){

  euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  ## Approximate thresholds
  # Find highest to low average metabolite
  avgMetabols = sort(apply(data_miss, 1, mean, na.rm=TRUE),
                     decreasing = TRUE, index = TRUE)
  # Sort vertically high to low metabolite
  rowSortedData = data_miss[avgMetabols$ix,]
  fullSortedData = split(rowSortedData, seq(nrow(rowSortedData)))

  fullSortedData = lapply(fullSortedData, function(x){
    sum(is.na(x)/ncol(data_miss))# Count missing vals per metabolite
  })
  fullSortedData = as.data.frame(do.call(rbind, fullSortedData))

  if (PercentMiss < 5){
    thresh_I = c(5)
  } else {
    thresh_I = seq(5,PercentMiss,5)
  }

  thresh_II = seq(60,80,5)
  thresh_III = seq(5, 60, 5)
  grid = base::expand.grid(thresh_I, thresh_II, thresh_III)
  threshs = list()
  distances = numeric()
  try({
    for (i in seq_len(nrow(grid))) {
      alpha = grid$Var1[i]
      beta = grid$Var2[i]
      gamma = grid$Var3[i]
      ## Approximate thresholds
      data_test = suppressWarnings({
        tryCatch({
          removeDataMM(data_sub, percentMiss = PercentMiss, alpha, beta, gamma)
        }, error = function(e) e)
      })

      if(inherits(data_test, "error")) next
      # Find highest to low average metabolite
      data_avgMetabols = sort(apply(data_test, 1, mean, na.rm=TRUE),
                              decreasing = TRUE, index = TRUE)
      # Sort vertically high to low metabolite
      data_rowSortedData = data_test[data_avgMetabols$ix,]
      data_fullSortedData = split(data_rowSortedData,
                                  seq(nrow(data_rowSortedData)))

      data_fullSortedData = lapply(data_fullSortedData, function(x){
        sum(is.na(x)/ncol(data_test))# Count missing vals per metabolite
      })
      data_fullSortedData = as.data.frame(do.call(rbind, data_fullSortedData))
      distance = euc.dist(data_fullSortedData, fullSortedData)
      threshs[[i]] = cbind(alpha, beta, gamma) # Store current threshold I
      distances[i] = distance # Store distance computed

    }
  }, silent = TRUE)

  return(list(thresh = threshs[which.min(distances)],
              distance = distances[which.min(distances)]))
}

# function to get largest complete data subset
largest_complete_subset = function(data_miss){
  shuffle_rows = function(x) {
    x = x[, sample(ncol(x))]
    return(x)
  }
  moveNAs = function(x) {
    # count NA
    num.na = sum(is.na(x))
    # remove NA
    x = x[!is.na(x)]
    # glue the number of NAs at the end
    x = c(x, rep(NA, num.na))
    return(x)
  }
  # Get largest possible subset of complete data
  data_original = shuffle_rows(data_miss)
  data_original = split(data_original, seq(nrow(data_miss)))
  data_original = lapply(data_original, moveNAs)

  SortedData = lapply(data_original, function(x){
    item = unlist(x)
    item = item[order(item)]
    return(item)

  })
  SortedData = do.call(rbind, SortedData)
  mostNAs = lapply(data_original, function(x){
    sum(is.na(x))
  })
  mostNAs.ix = which(unlist(mostNAs)==max(unlist(mostNAs)), arr.ind = TRUE)[1]
  data_original = do.call(rbind, data_original)
  remove.ix = ncol(data_original)-mostNAs[[mostNAs.ix]]
  data = data_original[,seq_len(remove.ix)]
  return(data)
}

# function to generate predictors/labels
generate_predictors = function(missingData, labels=TRUE){
  # Normalization function
  normaliz = function(x, eps=0.001){
    return((x-(min(x, na.rm = TRUE)-eps))/((max(x, na.rm = TRUE)-min(x, na.rm = TRUE))+eps))
    }

  numcols = ncol(missingData)
  numrows = nrow(missingData)
  # calc max_metabolite predictor
  max_metabolite = apply(apply(missingData, 2, as.numeric), 1, max, na.rm=TRUE)
  # calc min_metabolite predictor
  min_metabolite = apply(apply(missingData, 2, as.numeric), 1, min, na.rm=TRUE)
  # calc median_metabolite predictor
  median_metabolite = apply(apply(missingData, 2, as.numeric),
                            1, median, na.rm=TRUE)
  # calc average_metabolite predictor
  avg_metabolite = apply(apply(missingData, 2, as.numeric),
                         1, mean, na.rm=TRUE)
  # calc metabolite quantiles
  metabolite_quantiles = apply(apply(missingData, 2, as.numeric),
                               1, quantile, na.rm=TRUE)
  # calc amount missing per metabolite
  num_missing = apply(apply(missingData, 2, as.numeric),
                      1, function(x){sum(is.na(x))})
  # calc percent missing per metabolite predictor
  num_missing = num_missing/ncol(missingData)

  levels = matrix(NA, ncol = numcols, nrow = numrows) # initalize matrix
  # Use normalized data
  missingData_norm = t(apply(apply(missingData, 2, as.numeric), 1, normaliz))
  for (i in seq_len(numrows)){ # calc metabolite quantile categories predictor
    for (j in seq_len(numcols)){
      levels[i, j] = ifelse(missingData_norm[i,j] > metabolite_quantiles[3, i], "high",
                            ifelse(missingData_norm[i,j]<metabolite_quantiles[2, i], "low",
                                   ifelse(is.na(missingData_norm[i,j]), missingData[i,j],"medium")))
    }
  }
  levels[is.na(levels)] = "none"
  frame.as.vector = as.vector(t(missingData)) # Vectorize matrix

  # append max_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(max_metabolite, each=numcols))
  # append min_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(min_metabolite, each=numcols))
  # append median_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(median_metabolite, each=numcols))
  # append average_metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(avg_metabolite, each=numcols))
  # append percent missing per metabolite predictor
  frame.as.vector = cbind(frame.as.vector, rep(num_missing, each=numcols))
  # append metabolite quantile categories predictor
  frame.as.vector = cbind(frame.as.vector, as.vector(t(levels)))

  if (labels){
    target = ifelse(grepl("^[A-Za-z\\/]+$", frame.as.vector[,1], perl = TRUE),
                    frame.as.vector[,1], 'O') # Generate label vector
    frame.as.vector = cbind(frame.as.vector, target) # append labels
  }

  # normalize metabolite abundances
  missingvals = as.numeric(frame.as.vector[,1])
  # replace NAs with 0s
  missingvals[is.na(missingvals)] = 0
  # Append replacement to data frame
  frame.as.vector[,1] = missingvals
  # data structure change
  frame.as.vector = as.data.frame(frame.as.vector)
  # Ensure numeric predictors are numeric type
  frame.as.vector[,seq_len(6)] = apply(frame.as.vector[,seq_len(6)],
                                           2, as.numeric)
  # in case any rows are still missing omit
  frame.as.vector = na.omit(frame.as.vector)
  return(frame.as.vector)
}

# function to impute data usig specific imputation algorithms
imputation_algorithms = function(missingData,
                                 predictions,
                                 MCAR_algorithm,
                                 MNAR_algorithm,
                                 n_cores){

  # Imputation step
  numcols = ncol(missingData)
  numrows = nrow(missingData)
  predicted_missingData = matrix(predictions,
                                 ncol = numcols, nrow = numrows, byrow = TRUE)
  MCAR_index = which(predicted_missingData == "MCAR", arr.ind = TRUE)
  MNAR_index = which(predicted_missingData == "MNAR", arr.ind = TRUE)
  non_missing = which(predicted_missingData == "O", arr.ind = TRUE)
  Imputed_data = matrix(NA, ncol = numcols, nrow = numrows)

  if (n_cores != 1){
    plan(multisession, workers = 2)

  SingleAlgorithmImputations = future_lapply(list(MNAR_algorithm,
                                                  MCAR_algorithm),function(l){

    if (l == "Single"){
      MNAR_imputation = kapImpute2(missingData)
    }

    if (l == 'nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MNAR_imputation = nsKNN(missingData, k, iters=1)
    }

    if (l == "BPCA"){
      MCAR_imputation = pca(missingData,
                            method = "bpca",
                            verbose=FALSE)@completeObs
    }
    if (l == "random_forest"){
      MCAR_imputation = missForest(missingData)[["ximp"]]
    }
    if (l == 'Multi_nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MCAR_imputation = nsKNN(missingData, k,
                              iters=4, weighted=TRUE, scale = TRUE,
                              shuffle = TRUE)
    }

    if (exists('MNAR_imputation')){
      return(MNAR_imputation = MNAR_imputation)
    } else {
      return(MCAR_imputation = MCAR_imputation)
    }
  }, future.globals = c("kapImpute2", "nsKNN"), future.seed=TRUE)

  plan(sequential)


  MNAR_imputation = SingleAlgorithmImputations[[1]]
  MCAR_imputation = SingleAlgorithmImputations[[2]]
  } else {
    if (MNAR_algorithm == "Single"){
      MNAR_imputation = kapImpute2(missingData)
    }

    if (MNAR_algorithm == 'nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MNAR_imputation = nsKNN(missingData, k, iters=1)
    }

    if (MCAR_algorithm == "BPCA"){
      MCAR_imputation = pca(missingData,
                            method = "bpca",
                            verbose=FALSE)@completeObs
    }
    if (MCAR_algorithm == "random_forest"){
      MCAR_imputation = missForest(missingData)[["ximp"]]
    }
    if (MCAR_algorithm == 'Multi_nsKNN'){
      k = floor(sqrt(nrow(missingData)))
      if (k>numcols){
        k=numcols-1
      }
      MCAR_imputation = nsKNN(missingData, k,
                              iters=4, weighted=TRUE, scale = TRUE,
                              shuffle = TRUE)
    }
  }

  Imputed_data[MCAR_index] = MCAR_imputation[MCAR_index]
  Imputed_data[MNAR_index] = MNAR_imputation[MNAR_index]
  Imputed_data[is.na(Imputed_data)] = missingData[non_missing]
  return(list(MAI = Imputed_data,
              MCAR_imps = as.matrix(MCAR_imputation),
              MNAR_imps = as.matrix(MNAR_imputation)))
}
