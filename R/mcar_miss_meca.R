#'@title Impute NaNs in a simulated data according the MCAR mechanism
#'
#'@description Simulate NaNs in a complete data composed of discrete and continuous traits which are
#'correlated or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MCAR. It used the
#'function delete_MCAR from the package missMethods.
#'
#'@usage mcar_miss_meca(missingRate, ds, cols_mis)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds data frame in which NA should be created
#'@param cols_mis vector of index or columns name where missing values will be created
#'@return a dataset with NA following a pattern of MCAR.
#'@export
mcar_miss_meca <- function(missingRate, ds, cols_mis){

  if(is.numeric(cols_mis)){
    cols_mis <- names(ds)[cols_mis]
  }

  #call MCAR function
  if(ncol(ds) == 1){
    missingMCAR <- missMethods::delete_MCAR(ds, missingRate, cols_mis = cols_mis)
  }

  if(ncol(ds) > 1){
    missingMCAR <- missMethods::delete_MCAR(ds, missingRate, cols_mis = cols_mis, p_overall = FALSE)
  }

  #check if all the states are represented. (discrete traits)
  #discIndex <- grep("I.", names(ds[, cols_mis, drop = FALSE]))
  naColumns <- which(colSums(is.na(ds[, cols_mis, drop = FALSE])) > 0)
  classVar <- lapply(ds[, cols_mis, drop = FALSE], class)
  discIndex <- which(classVar != "numeric")


  if(length(discIndex) != 0){

    for(t in 1:length(discIndex)){

      nbrStates <- unique(ds[, discIndex[t]])
      nbrStatesMiss <- unique(na.omit(missingMCAR[, discIndex[t]]))
      diff <- setdiff(nbrStates, nbrStatesMiss)

      #if > 0 means not all the states are represented.
      if(length(diff) != 0){

        rowIndexSave <- c()

        for(s in 1:length(diff)){

          rowIndexSate <- which(ds[, discIndex[t]] == diff[s])

          indexToKeep <- sample(rowIndexSate, 1)

          rowIndexSave <- c(rowIndexSave, indexToKeep)

        }

        newMissingRate <- missingRate * nrow(ds[-rowIndexSave, ]) / nrow(ds)

        #call MCAR function
        newMissingMCAR <- missMethods::delete_MCAR(ds[-rowIndexSave, ], newMissingRate, cols_mis = cols_mis, p_overall = FALSE)

        missingMCAR <- ds
        missingMCAR[rownames(newMissingMCAR), ] <- newMissingMCAR
      }
    }
  }
  return(missingMCAR)
}
