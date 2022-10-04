#' @title non-parametric missing value imputation for mixed-type data by kNN(k-nearest neighbor algorithm)
#'
#' @description This function imputes missing data for continuous and discrete traits applying the kNN approach
#'
#' @usage kNN_phylo(missingData, k, numFun, catFun, variance_fraction = 0, tree, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param k integer, number of nearest neighbours used
#' @param numFun numFun: function for aggregating the kNN in the case of a numerical variable
#' @param catFun catFun: function for aggregating the kNN in the case of a categorical variable
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param tree phylo object
#' @param hint data frame already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return A data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' @export
kNN_phylo <- function(missingData, k, numFun, catFun, variance_fraction = 0, tree, hint = NULL){

  NbrCol <- ncol(missingData)

  #include imputed data in case
  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))

    missingData <- cbind(missingData, hint)
  }

  if(variance_fraction != 0 & variance_fraction != 2){

    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }


  #get the column name with NaN
  variable <- c()
  for(c in 1:ncol(missingData)){
    if(any(is.na(missingData[,c]))){
      variable <- c(variable, colnames(missingData)[c])
    }
  }

  DataImputed <- VIM::kNN(missingData, k = k, variable = variable,
                          numFun = numFun, catFun = catFun, imp_var = FALSE)[,1:NbrCol,  drop = FALSE]
  #keep the tip names
  rownames(DataImputed) <- rownames(missingData)

  parameters <- list(k = k, numFun = numFun, catFun = catFun)

  return(list(imputedData = DataImputed, parameters = parameters))
}
