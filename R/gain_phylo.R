#' @title GAIN
#'
#' @description This function imputes missing data for continuous and discrete traits applying the Generative Adversarial #'
#' Imputation Networks (GAIN)
#'
#' @usage gain_phylo(missingData, variance_fraction, Data, batch_size = round(ncol(missingData)*0.2), hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL)
#'
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the
#' phylogenetic inertia
#' @param tree phylo object
#' @param batch_size integer
#' @param hint_rate numerical
#' @param alpha numerical, hyperparameter
#' @param epochs integer, iterations
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return A list of list containing in the "tab" imputedData the imputed Data and in the "tab" parametersGain, the discriminant loss values, the generative loss values, the MSE loss values and the iterations correspond to these values (important to generate a plot).
#' @export
gain_phylo <- function(missingData, variance_fraction, tree, batch_size = round(ncol(missingData)*0.2),
                  hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL){


  #load python files
  reticulate::source_python(system.file("python", "modelV2R.py", package = "TDIP"))

  NbrCol <- ncol(missingData)

  #get the factors columns
  factorsColumns <- names(Filter(is.factor, missingData))

  rNames <- row.names(missingData)
  colNames <- colnames(missingData)

  #include imputed data in case
  if(!is.null(hint)){
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

  #convert factors as dummy
  if(length(factorsColumns) != 0){
    oneHotConvertion <- generateOneHotEncoVariables(missingData)
    missingData <- oneHotConvertion$data
    oneHotColnames <- colnames(missingData)
  }

  #convert all the values as numeric (second check)
  missingData <- apply(missingData, 2, as.numeric)

  #run Gain
  obj <- Gain()
  gainOutput <- obj$runGain(missingData, batch_size, hint_rate, alpha, epochs)

  #merge list of list
  for (l in 2:length(gainOutput)){
    gainOutput[[l]] <- unlist(gainOutput[[l]], recursive = FALSE)
  }

  #convert imputedData in dataframe
  imputedData <- as.data.frame(gainOutput[[1]])

  if(length(factorsColumns) != 0){
    names(imputedData) <- oneHotColnames
    imputedData <- convertOneHotEncoIntoFactors(imputedData, colNames)
  }

  #remove hint variables or phylogenetic eigenvectors (in case no one hot encoding)
  imputedData <- imputedData[, 1:NbrCol, drop = FALSE]
  names(imputedData) <- colNames
  row.names(imputedData) <- rNames

  #parameters
  lossValues <- gainOutput[c(2:5)]
  names(lossValues) <- c("d_loss", "g_loss", "mse_loss", "epochs")

  return(list(imputedData = imputedData, parameters = lossValues))
}
