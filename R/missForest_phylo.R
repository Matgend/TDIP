#' @title Non-parametric missing values imputation for mixed-type data by missForest(randomForest)
#'
#' @description This function imputes missing data for continuous and discrete traits applying the missForest approach
#'
#' @usage missForest_phylo(missingData, variance_fraction = 0, maxiter = 10, ntree = 100,
#' mtry = sqrt(ncol(missingData)), Data, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param maxiter maximum number of iterations to be performed given the stopping criterion is not met beforehand.
#' @param ntree number of trees to grow in each forest.
#' @param mtry number of variables randomly sampled at each split. By default it's the square root of the number of
#' variables
#' @param tree phylo object
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return A data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' @export
#'
missForest_phylo <- function(missingData, variance_fraction = 0, maxiter = 10, ntree = 100,
                             mtry = sqrt(ncol(missingData)), tree, hint = NULL){

  colNames <- names(missingData)


  #missForest doesn't support when more than 53 classes then remove the columns having (is no NA inside)
  splitData <- NULL
  for(col in ncol(missingData):1){
    if(is.factor(missingData[,col])){
      if(length(levels(missingData[,col])) > 52){
        if(any(is.na(missingData[,col]))){
          stop("missForest doesn't able to impute missing value in a variable with more than 53 categories")
        }
        else{

          if(is.null(splitData)){
            splitData <- missingData[,col, drop = FALSE]
          }

          else{
            splitData <- cbind(splitData, missingData[,col, drop = FALSE])
          }

          missingData[,col] <- NULL

        }
      }
    }
  }

  Nvariables <- ncol(missingData)

  #include imputed data in case
  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))

    missingData <- cbind(missingData, hint)
  }

  # want to include phylogeny information
  if(variance_fraction != 0 & variance_fraction != 2){

    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }

  #run missForest
  #missForest_imputation <- missForest::missForest(xmis = missingData, maxiter = maxiter,
  #                                                ntree = ntree, mtry = mtry, verbose = TRUE)

  missForest_imputation <- myMissForest(xmis = missingData, maxiter = maxiter,
                                        ntree = ntree, mtry = mtry, verbose = TRUE)

  myMissForest
  #cut eigenvectors columns
  missForest_imputation$ximp <- missForest_imputation$ximp[,1:Nvariables, drop = FALSE]


  if(!is.null(splitData)){
    missForest_imputation$ximp <- cbind(missForest_imputation$ximp, splitData)
    missForest_imputation$ximp <- missForest_imputation$ximp[, colNames]
  }

  parameters <- list(maxiter = maxiter, ntree = ntree, mtry = mtry)

  return(list(imputedData = missForest_imputation$ximp, parameters = parameters))
}
