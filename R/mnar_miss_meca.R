#'@title Impute NaNs in a simulated data according the MNAR mechanism
#'
#'@description Simulate NaNs in a complete data composed of discrete and continuous traits which are correlated
#'or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MNAR. It used the function
#'delete_MNAR from the package missMethods.
#'
#'@usage mnar_miss_meca(missingRate, ds, cols_mis)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds dataframe in which NA should be created
#'@param cols_mis vector of index or columns name where missing value will be created
#'@return a dataset with NA following a pattern of MNAR.
#'@export
mnar_miss_meca <- function(missingRate, ds, cols_mis){

  if(is.numeric(cols_mis)){
    colPartition <- cols_mis
  }

  else if(is.character(cols_mis)){
    colPartition <- match(cols_mis, names(ds))
  }

  for(col in colPartition){

    naColumns <- which(colSums(is.na(ds[, col, drop = FALSE])) > 0)
    classVar <- lapply(ds[, col, drop = FALSE], class)
    contiIndex <- which(classVar == "numeric")
    if(length(contiIndex) == 1){
      ds[ , col] <- missMethods::delete_MNAR_censoring(ds[ , col, drop = FALSE], missingRate,
                                                       cols_mis = cols_mis[col], where = "upper") #don't change where arg
    }

    else{
      levelsVector <- sort(unique(ds[ ,col]))

      nbrValueToRemove <- round(nrow(ds) * missingRate)

      maxState <- tail(levelsVector, 1)

      #in case not enough value in max state, remove in the other states
      while (nbrValueToRemove > 0){

        levelsVector <- levelsVector[-which(levelsVector == maxState)]
        indexLargeState <- which(ds[,col] == maxState)

        #keep present 1 value of the state
        indexLargeState <- indexLargeState[-which(indexLargeState == sample(indexLargeState, 1))]

        if(length(indexLargeState) >= nbrValueToRemove){
          indexLargeState <- sample(indexLargeState, nbrValueToRemove)
          nbrValueToRemove <- 0

        }

        else{
          nbrValueToRemove <- nbrValueToRemove - length(indexLargeState)
          indexLargeState <- sample(indexLargeState, length(indexLargeState))
          maxState <- sample(levelsVector, 1)
        }

        #apply NA
        ds[indexLargeState, col] <- NA
      }
    }
  }
  return(ds)
}
