#' @title Imputation of missing values in mixed data by the MICE approach
#'
#' @description This function imputes missing data for continuous and discrete traits applying the MICE approach using the
#' predictive mean matching method for continuous traits, the polytomous regression for unordered categorical traits having
#' more than 2 states, if binary the algorithm uses a logistic regression.
#'
#' @usage mice_phylo(missingData, nbrMI, variance_fraction = 0, tree, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param nbrMI integer, mentioning the total number of imputations
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param tree phylo object
#' @param hint data frame already imputed by a comparative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' @export
mice_phylo <- function(missingData, nbrMI, variance_fraction = 0, tree, hint = NULL){

  colNames <- names(missingData)
  nativeMissingData <- missingData

  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))

    missingData <- cbind(missingData, hint)
  }

  if(variance_fraction != 0 & variance_fraction != 2){

    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[, 1:ncol(eigen), drop = FALSE])

  }
  names(missingData) <- paste0("A",as.character(1:ncol(missingData)))

  #in case collinearity is too important, and MICE can't impute the data
  ImputedMICE <- tryCatch(
    {

      #check the number of levels for factors
      levelsByColumns <- lapply(missingData, levels)
      lengthLevels <- lapply(levelsByColumns, length)
      columnsCatPMM <- as.numeric(lengthLevels > 20) #in literature said that in case nbr state > 20 use pmm instead of polyreg
      columnsLogReg <- as.numeric(lengthLevels == 2)

      method <- NULL

      if(sum(columnsCatPMM) != 0){
        method <- rep("polyreg", length(columnsCatPMM))
        method[which(columnsCatPMM == 1 | as.numeric(lengthLevels == 0) == 1)] <- "pmm"
        method[which(columnsLogReg == 1)] <- "logreg"
      }

      # ImputedMICE <- mice::mice(missingData, m = nbrMI, maxit = 5, method = method,
      #                           printFlag = FALSE,
      #                           pred = mice::quickpred(missingData, mincor = 0.1, minpuc = 0, include = "", exclude = "",
      #                                            method = "pearson"))

      ImputedMICE <- mice::mice(missingData, m = nbrMI, maxit = 5, method = method,
                                printFlag = FALSE)


      #choose the first column
      imputedData <- mice::complete(ImputedMICE, action = 1L)[, 1:length(colNames),  drop = FALSE]
      names(imputedData) <- colNames

      colWithNaN <-  which(colSums(is.na(imputedData)) != 0)
      if(length(colWithNaN) != 0){
        for(v in names(colWithNaN)){
          ry <- !is.na(imputedData[ ,v])
          print(paste("random imputation of column:", v))
          imputedData[which(!ry) ,v] <- mice::mice.impute.sample(imputedData[ ,v], ry)
        }
      }

      parameters <- list(nbrMI = nbrMI)
      list(imputedData = imputedData, parameters = parameters)
    },
    error = function(cond){
      print("MICE doesn't work")
      list(imputedData = nativeMissingData, parameters = "Collinearity induces imputation error")
    }
  )
  return(ImputedMICE)
}
