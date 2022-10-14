#' @title Calculate error of imputation
#'
#' @description This function calculates the RMSE for imputed continuous data, the absolute
#' error for imputed ordinal data and the Proportion of falsely classified entries (PFC) for the other
#' subcategories of discrete data
#'
#' @usage imputation_error(imputedData, trueData, missingData, imputationApproachesName, dataframe)
#'
#' @param imputedData array of imputed data
#' @param trueData array of true data
#' @param missingData array of data with missing values
#' @param imputationApproachesName character providing the name of the imputed appraoch used
#' @param dataframe dataframe given information about the data (type of variables, model of evolution, ....). If NULL, will measure the error of ordinate variables as categorical variables.
#' @return A data frame with in the first column the trait names and in the second the errors
#' @export
imputation_error <- function(imputedData, trueData, missingData, imputationApproachesName, dataframe = NULL){

  #get the ordinal trait reference
  ordinalTraits <- NULL
  if(!is.null(dataframe))
    ordinalTraits <- which(dataframe$class == "ordinal") #give the row in dataframe which correspond to /n in data names
  errors <- c()
  traitNames <- c()
  for (c in 1:ncol(missingData)){

    #know is NaNs in the columns(trait)
    NaNRowIndex <- which(is.na(missingData[,c]))

    if(length(NaNRowIndex != 0)){

      traitNames <- c(traitNames, names(trueData)[c])
      #missingValues <- missingData[NaNRowIndex, c]
      trueValues <- trueData[NaNRowIndex, c]
      imputedValues <- imputedData[NaNRowIndex, c]

      #in case continuous data
      classVar <- class(trueData[,c])

      contiColumns <- which(classVar == "numeric")

      if(length(contiColumns) != 0){
        #rmse
        error <- sqrt(mean((as.numeric(imputedValues) - as.numeric(trueValues))^2))

      }

      #in case ordinal trait
      else if(length(ordinalTraits) != 0 & length(grep(paste0("/", ordinalTraits), names(missingData)[c])) == 1){
        #imputation error for ordinal traits (absolute error)
        error <- mean(abs((as.numeric(imputedValues) - as.numeric(trueValues)) / as.numeric(trueValues)))
      }

      #in case discrete data
      else{
        error <- (sum(as.character(imputedValues) != as.character(trueValues)) / length(trueValues))
        #error <- length(setdiff(imputedValues, trueValues)) / length(trueValues)
      }
      errors <- c(errors, error)
    }
  }
  output <- data.frame(trait = traitNames, c2 = errors)
  names(output)[2] <- imputationApproachesName

  return(output)
}
