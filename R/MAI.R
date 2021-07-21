MAI = function(data_miss,
               MCAR_algorithm = "BPCA",
               MNAR_algorithm = "Single",
               n_cores = 1){

  # Possible imputation algorithms
  MCAR_algorithms = c("BPCA", "random_forest", "Multi_nsKNN")
  MNAR_algorithms = c("nsKNN", "Single")

  # Check argument inputs
  if (!MCAR_algorithm %in% MCAR_algorithms){
    stop("Incorrect MCAR algorithm selected. Possible choices include:BPCA,
         random_forest, or Multi_nsKNN ")
  }

  if (!MNAR_algorithm %in% MNAR_algorithms){
    stop("Incorrect MNAR algorithm selected. Possible choices include: nsKNN
         or Single ")
  }

  # Check data
  if (is(data_miss, "SummarizedExperiment")){
    if (sum(is.na(assay(data_miss))) == 0){
      stop("No missing values detected to impute")
    }
    if (any(rowSums(is.na(assay(data_miss))) == ncol(data_miss))){
      stop("Detected a row that contains no observed data, omit this row.")
    }
    if (any(rowSums(is.na(assay(data_miss))) > 0.8*ncol(data_miss))){
      warning("Detected a row with more than 80% missing values.
            Consider omitting rows that are > 80% missing.")
    }
  }else{

  if (sum(is.na(data_miss)) == 0){
    stop("No missing values detected to impute")
  }
  if (any(rowSums(is.na(data_miss)) == ncol(data_miss))){
    stop("Detected a row that contains no observed data, omit this row.")
  }
  if (any(rowSums(is.na(data_miss)) > 0.8*ncol(data_miss))){
    warning("Detected a row with more than 80% missing values.
            Consider omitting rows that are > 80% missing.")
  }}

  # Check n_cores
  if (n_cores%%1 != 0){
    stop("n_cores must be an integer greater than or equal to -1.")
  }
  if (n_cores > detectCores()){
    stop("You specified more cores than is available to your system.")
  }
  if (n_cores < -1){
    stop("n_cores must be an integer greater than or equal to -1.")
  }

  if (is(data_miss, "SummarizedExperiment")) {
    fldata = data_miss
    data_miss = assay(data_miss)
  }

  # Original data percent missing (total)
  PercentMiss = (sum(is.na(data_miss))/(nrow(data_miss)*ncol(data_miss)))*100
  data = largest_complete_subset(data_miss)
  # Estimate threshold I
  threshs = list()
  print("Estimating pattern of missingness")
  for (i in seq_along(1:10)){
    try({threshs[[i]] = check_distance(data_miss, data, PercentMiss)},
        silent = FALSE)
  }
  smallestDistance = threshs[[which.min(lapply(threshs, function(x){
    unlist(x[["distance"]])
  }))]][["thresh"]][[1]]
  # Set Threshold params
  alpha = smallestDistance[1]
  beta = smallestDistance[2]
  gamma = smallestDistance[3]


  # generate missingness
  print('Imposing missingness')
  missingData = removeDataMM_altered(data,
                                     percentMiss = PercentMiss,
                                     percentBelowThresh_I = alpha,
                                     LowAbundThresh_II = beta,
                                     percentMisslowAbund_III = gamma)

  # generate predictors
  print('Generating features')
  full_data = suppressWarnings(generate_predictors(missingData, labels=TRUE))

  # Split into training and testing sets
  trainIndex = createDataPartition(full_data$target,
                                   p=0.95,
                                   list = FALSE,
                                   times = 1)
  trainingSet = full_data[trainIndex,]
  testingSet = full_data[-trainIndex,]
  # Train model
  fitControl = trainControl(
    method = "cv",
    number = 2)
  # Parallelize training
  if (n_cores > 1)
  {
    cores = n_cores
    cl = makePSOCKcluster(cores)
    registerDoParallel(cl)
  }
  if (n_cores == -1)
  {
    cores = detectCores()
    cl = makePSOCKcluster(cores)
    registerDoParallel(cl)
  }

  print("Training")
  # Random forest
  rfFit = train(target ~ .,
                data = trainingSet,
                method = "rf",
                trControl = fitControl,
                ntree = 300,
                verbose = FALSE)

  if (n_cores != 1){
    stopCluster(cl)
  }


  # Use model to make predictions on initial data
  print("Predicting")
  predictions = predict(rfFit,
                        newdata = suppressWarnings(
                          generate_predictors(data_miss, labels=FALSE))
                        )
  # Impute based on classification predictions
  print("Imputing")
  Imputed_data = imputation_algorithms(data_miss,
                                       predictions,
                                       MCAR_algorithm = MCAR_algorithm,
                                       MNAR_algorithm = MNAR_algorithm,
                                       n_cores = n_cores)

  if(exists("fldata")){
    assay(fldata, 1, withDimnames = FALSE) = as.matrix(Imputed_data[[1]])
    MCAR_imputations = Imputed_data[[2]]
    MNAR_imputations = Imputed_data[[3]]
    estimated_Params = list(
      Alpha = alpha,
      Beta = beta,
      Gamma = gamma
    )
    metadata(fldata) = list(metadata(fldata),
                            Estimated_Params = estimated_Params,
                            MCAR_imputations = MCAR_imputations,
                            MNAR_imputations = MNAR_imputations)
    return(fldata)
  }else{
    return(list(Imputed_data = Imputed_data,
                Estimated_Params = list(
                  Alpha = alpha,
                  Beta = beta,
                  Gamma = gamma
                )))
  }


}
