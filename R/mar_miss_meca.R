#'@title Impute NaNs in a simulated data according the MAR mechanism
#'
#'@description Simulate NaNs in a complete data composed of discrete and continuous traits which are correlated
#'or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MAR. It used the function delete_MAR
#'from the package missMethods.
#'
#'@usage mar_miss_meca(missingRate, ds, cols_mis, cols_ctrl)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds dataframe in which NA should be created
#'@param cols_mis vector of index or columns name where missing value will be created
#'@param cols_ctrl vector of index or columns name which control the creation of missing values (same size of cols_mis)
#'@return A dataset with NA following a pattern of MAR.
#'@export
mar_miss_meca <- function(missingRate, ds, cols_mis, cols_ctrl){

  if(is.numeric(cols_mis)){
    cols_mis <- names(ds)[cols_mis]
  }

  #call MAR function
  missingMAR <- missMethods::delete_MAR_censoring(ds, missingRate,
                                                  cols_mis = cols_mis, cols_ctrl = cols_ctrl,
                                                  where = "upper")

  #check if all the states are represented. (discrete traits)

  #discIndex <- grep("I.", names(ds[, cols_mis, drop = FALSE]))
  naColumns <- which(colSums(is.na(ds[, cols_mis, drop = FALSE])) > 0)
  classVar <- lapply(ds[, cols_mis, drop = FALSE], class)
  discIndex <- which(classVar == "factor")

  if(length(discIndex) != 0){

    for(t in 1:length(discIndex)){
      nbrStates <- unique(ds[, discIndex[t]])

      nbrStatesMiss <- unique(na.omit(missingMAR[, discIndex[t]]))

      diff <- setdiff(nbrStates, nbrStatesMiss)

      #if > 0 means not all the states are represented.
      if(length(diff) != 0){

        rowIndexSave <- c()

        for(s in 1:length(diff)){

          #index missing States
          indexMissStates <- which(ds[, discIndex[t]] == diff[s])

          #isolate values of ctrl
          valuesCor <- ds[, cols_ctrl[t]][indexMissStates]

          #if continuous values
          if(is.numeric(valuesCor)){

            rowIndexMinVal <- which.min(valuesCor)
            rowIndexSave <- c(rowIndexSave, indexMissStates[rowIndexMinVal])
          }

          else{
            rowIndexMinVal <- sample(indexMissStates, 1)
            rowIndexSave <- c(rowIndexSave, rowIndexMinVal)
          }
        }

        newMissingRate <- missingRate * nrow(ds[-rowIndexSave, ]) / nrow(ds)
        #call MAR function
        newMissingMAR <- missMethods::delete_MAR_censoring(ds[-rowIndexSave, ], newMissingRate,
                                                           cols_mis = cols_mis, cols_ctrl = cols_ctrl,
                                                           where = "upper")

        missingMAR <- ds
        missingMAR[rownames(newMissingMAR), ] <- newMissingMAR
      }
    }
  }
  return(missingMAR)
}
