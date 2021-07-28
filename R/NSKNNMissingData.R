removeDataMM = function(data,
                        percentMiss,
                        percentBelowThresh_I,
                        LowAbundThresh_II,
                        percentMisslowAbund_III){
  data = as.matrix(data)
  percentMisslowAbund_III = percentMisslowAbund_III/0.9
  percentMVmedAbund = 0.5*percentMisslowAbund_III
  TotalNumValRemove = round(percentMiss/100*ncol(data)*nrow(data))
  data[data == 0] = 0.00001

  temp = sort(apply(data, 1, mean), decreasing = FALSE, index = TRUE)
  avgMetabolite = temp[["x"]]
  sortAvgIdx = temp[["ix"]]
  numMetBelowThresh = round(percentBelowThresh_I/100*length(avgMetabolite))
  belowThreshMet = data[sortAvgIdx[seq_len(numMetBelowThresh)],]
  belowThreshIdx = sortAvgIdx[seq_len(numMetBelowThresh)]
  aboveThreshMet = data[sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)],]
  aboveThreshIdx = sortAvgIdx[(numMetBelowThresh+1):length(avgMetabolite)]

  numMetLowAbund = round(LowAbundThresh_II/100*nrow(belowThreshMet))
  lowAbundMet = as.matrix(data[sortAvgIdx[seq_len(numMetLowAbund)],])
  medAbundmet = as.matrix(
    data[sortAvgIdx[(numMetLowAbund+1):nrow(belowThreshMet)],]
    )
  # # of samples for each metabolite that belong to the low abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshLowAbund = round(
    percentMisslowAbund_III/100*ncol(lowAbundMet)
    )
  if (numSampBelowAbundThreshLowAbund == 0){dw = 0}else{dw=1}
  # Sort matrix with low mean abundance metabolites from low to high abundance
  # values for each metabolite.
  sortLowAbundmet = t(apply(lowAbundMet, 1, sort))
  sortAbundLowAbundidx = t(apply(lowAbundMet, 1, order))

  # Create matrix with low abundance values for each metabolite with a low mean
  # abundance.
  LowLowAbund = as.matrix(
    sortLowAbundmet[,seq_len(numSampBelowAbundThreshLowAbund)]
    )
  # Create a matrix with high abundance values for each metabolite with a low
  # mean abundance
  HighLowAbund = as.matrix(
    sortLowAbundmet[,(numSampBelowAbundThreshLowAbund+1):ncol(sortLowAbundmet)]
    )

  alwaysMV = round(0.8*ncol(LowLowAbund))
  LowLowAbund[,seq_len(alwaysMV)] = NA

  if (alwaysMV < ncol(LowLowAbund)){
    LowLowAbundRandomMV = as.matrix(
      LowLowAbund[,(alwaysMV+1):ncol(LowLowAbund)]
      )
    LowLowAbund = as.matrix(LowLowAbund[,-((alwaysMV+1):ncol(LowLowAbund))])
    numVal = ncol(LowLowAbundRandomMV)*nrow(LowLowAbundRandomMV)
    numValRemove = round(0.5*numVal)
    removeIdx = sample(numVal,numValRemove)
    temp = as.vector(t(LowLowAbundRandomMV))
    temp[removeIdx] = NA
    LowLowAbundRandomMV = matrix(temp, nrow = nrow(LowLowAbundRandomMV),
                                 ncol = ncol(LowLowAbundRandomMV), byrow = TRUE)
    LowLowAbund = cbind(LowLowAbund, LowLowAbundRandomMV)
  }

  # # of samples for each metabolite that belong to the medium abundance group
  # of metabolites that will be missing.
  numSampBelowAbundThreshMedAbund = round(
    percentMVmedAbund/100*ncol(medAbundmet)
    )
  if (numSampBelowAbundThreshMedAbund == 0){dz = 0}else{dz=1}
  # Sort matrix with high mean abundance metabolites from low to high abundance
  # values for each metabolite.
  sortmedAbundmet = t(apply(medAbundmet,1, sort))
  sortAbundMedAbundidx = t(apply(medAbundmet,1, order))
  # Create matrix with low abundance values for each metabolite with a high mean
  # abundance.
  LowMedAbund = as.matrix(
    sortmedAbundmet[,seq_len(numSampBelowAbundThreshMedAbund)]
    )
  # Create a matrix with high abundance values for each metabolite with a high
  # mean abundance.
  HighMedAbund = as.matrix(
    sortmedAbundmet[,(numSampBelowAbundThreshMedAbund+1):ncol(sortmedAbundmet)]
    )
  alwaysMV = round(0.8*ncol(LowMedAbund))
  LowMedAbund[,seq_len(alwaysMV)] = NA

  if (alwaysMV < ncol(LowMedAbund)){
    LowMedAbundRandomMV = as.matrix(
      LowMedAbund[,(alwaysMV+1):ncol(LowMedAbund)]
      )
    LowMedAbund = as.matrix(LowMedAbund[,-((alwaysMV+1):ncol(LowMedAbund))])
    numVal = ncol(LowMedAbundRandomMV)*nrow(LowMedAbundRandomMV)
    numValRemove = round(0.5*numVal)
    removeIdx = sample(numVal,numValRemove)
    temp = as.vector(t(LowMedAbundRandomMV))
    temp[removeIdx] = NA
    LowMedAbundRandomMV = matrix(temp,
                                 nrow = nrow(LowMedAbundRandomMV),
                                 ncol = ncol(LowMedAbundRandomMV),
                                 byrow = TRUE)
    LowMedAbund = cbind(LowMedAbund, LowMedAbundRandomMV)
  }

  # Add MCAR values
  numMVbelowThreshLOD = (sum(sum(is.na(LowMedAbund))) +
                           sum(sum(is.na(LowLowAbund))))*dz
  numMVrandom = TotalNumValRemove - numMVbelowThreshLOD
  if (numMVrandom < 0){
    numMVrandom = 0
  }

  MVrandomAvailable = ncol(data)*nrow(data)-dz*ncol(LowMedAbund)*nrow(LowMedAbund)-dw*ncol(LowLowAbund)*nrow(LowLowAbund)
  aboveThreshWithMV = aboveThreshMet
  HighLowAbundWithMV = HighLowAbund
  HighMedAbundWithMV = HighMedAbund
  size_aboveThreshWithMV = dim(aboveThreshWithMV)
  size_HighLowAbundWithMV = dim(HighLowAbundWithMV)
  size_HighMedAbundWithMV = dim(HighMedAbundWithMV)

  remainingValuesVector = c(as.vector(t(aboveThreshWithMV)),
                            as.vector(t(HighLowAbundWithMV)),
                            as.vector(t(HighMedAbundWithMV)))
  removeIdx = sample(length(remainingValuesVector), numMVrandom)
  remainingValuesVector[removeIdx] = NA

  aboveThreshWithMV = remainingValuesVector[seq_len((ncol(aboveThreshWithMV)*nrow(aboveThreshWithMV)))]
  HighLowAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+1):(length(aboveThreshWithMV)+ncol(HighLowAbundWithMV)*nrow(HighLowAbundWithMV))]
  HighMedAbundWithMV = remainingValuesVector[(length(aboveThreshWithMV)+(length(HighLowAbundWithMV)+1)):length(remainingValuesVector)]

  aboveThreshWithMV = matrix(aboveThreshWithMV,
                             nrow = size_aboveThreshWithMV[1],
                             ncol = size_aboveThreshWithMV[2],
                             byrow = TRUE)
  HighLowAbundWithMV = matrix(HighLowAbundWithMV,
                              nrow = size_HighLowAbundWithMV[1],
                              ncol = size_HighLowAbundWithMV[2],
                              byrow = TRUE)
  HighMedAbundWithMV = matrix(HighMedAbundWithMV,
                              nrow = size_HighMedAbundWithMV[1],
                              ncol = size_HighMedAbundWithMV[2],
                              byrow = TRUE)

  # Recreate original matrix with missing values and use sortrows to
  # reorder values to their original indices.
  if(dz==0){
    MedAbundWithMV = HighMedAbundWithMV
  }else{
    MedAbundWithMV = cbind(LowMedAbund, HighMedAbundWithMV)
  }
  #MedAbundWithMV = apply(MedAbundWithMV, 2, sort, na.last=T)
  if (ncol(MedAbundWithMV) > 1){
    for (i in seq_len(nrow(MedAbundWithMV))){
      MedAbundWithMV[i,] = MedAbundWithMV[i,][order(sortAbundMedAbundidx[i,])]
    }}

  if(dw==0){
    LowAbundWithMV = HighLowAbundWithMV
  }else{
    LowAbundWithMV = cbind(LowLowAbund, HighLowAbundWithMV)
  }
  if (ncol(LowAbundWithMV)>1){
    for (i in seq_len(nrow(LowAbundWithMV))){
      LowAbundWithMV[i,] = LowAbundWithMV[i,][order(sortAbundLowAbundidx[i,])]
    }}
  if (dz==0 & dw == 0){
    belowThreshWithMV = cbind(LowAbundWithMV, MedAbundWithMV)
  } else if (dz==0 & dw ==1){
    belowThreshWithMV = rbind(LowAbundWithMV, t(MedAbundWithMV))
  }else{
    belowThreshWithMV = rbind(LowAbundWithMV, MedAbundWithMV)
  }
  if (dz==0 & dw==0){
    dataMVtemp = rbind(t(belowThreshWithMV), aboveThreshWithMV)
  }else{
    dataMVtemp = rbind(belowThreshWithMV, aboveThreshWithMV)
  }

  dataMV = dataMVtemp[order(sortAvgIdx),]
  return(dataMV)
}
