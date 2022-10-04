#' @title Impute dataframe containing missing values
#'
#' @description This function imputes the dataset containing missing values with Rphylopars, corHMM,
#' MICE, missForest, KNN and GAIN according to 3 strategies, without phylogeny, with phylogenetic informations and 2-step
#'
#' @usage missing_data_imputation(ImputationApproachesNames, data, tree = NULL, strategies, varfrac = NULL, save = NULL)
#'
#' @param ImputationApproachesNames names of the imputations methods (character vector)
#' @param data dataframe containing NAs
#' @param tree phylo object
#' @param strategies character vector "NP", "P", "2-step"
#' @param varfrac amount of variance explained by the eigenvectors of the phylogenetic tree (provide only if in strategy
#' there is "P")
#' @param save name of the output
#' @return return a .RData object with the data imputed and all the parameters used for.
#' @export
missing_data_imputation <- function(ImputationApproachesNames, data, tree = NULL, strategies,
                            varfrac = NULL, save = NULL){

  #variables to define for imputation
  nbrMI <- 1
  k <- 2
  numFun <- laeken::weightedMedian
  catFun <- VIM::maxCat
  hint <- NULL

  strategies <- sort(strategies) #order is "2-step", "NP", "P"

  imputedData <- vector()
  imputedNames <- c()

  naColumns <- which(colSums(is.na(data)) > 0)
  classVar <- lapply(data, class)
  contiColumns <- which(classVar == "numeric")
  discColumns <- which(classVar != "numeric")
  colNames <- names(data)


  #PI
  if(("pi_continuous_traits" %in% ImputationApproachesNames | "pi_categorical_traits" %in% ImputationApproachesNames &
     "P" %in% strategies) | "2-step" %in% strategies){

    imputedValue <- data

    contiImputed <- NULL
    discrImputed <- NULL

    if(any(contiColumns %in% naColumns) & ("pi_continuous_traits" %in% ImputationApproachesNames | "2-step" %in% strategies)){
      contiImputed <- pi_continuous_traits(data[ ,contiColumns[which(contiColumns %in% naColumns)], drop = FALSE], tree)
      #contiImputed$imputedData <- exp(contiImputed$imputedData)
      imputedValue[ ,contiColumns] <- contiImputed$imputedData
      print("done Rphylopars")
    }



    if(any(discColumns %in% naColumns) & ("pi_categorical_traits" %in% ImputationApproachesNames | "2-step" %in% strategies)){
      discrImputed <- pi_categorical_traits(data[ ,discColumns[which(discColumns %in% naColumns)], drop = FALSE], tree)
      imputedValue[, discColumns] <- discrImputed$imputedData
      print("done corHMM")
    }


    if("P" %in% strategies){
      names(imputedValue) <- colNames
      imputedData <- c(imputedData, list(imputedValue))
      imputedName <- "PI"
      imputedNames <- c(imputedNames, imputedName)
    }

    if("2-step" %in% strategies){
      if(!is.null(contiImputed) & !is.null(discrImputed)){
        hint <- cbind(contiImputed$imputedData, discrImputed$probabilities)
      }
      else if(!is.null(contiImputed) & is.null(discrImputed)){
        hint <- contiImputed$imputedData
      }
      else{
        hint <- discrImputed$probabilities
      }

    }

  }

  print("start ML")
  for(s in 1:length(strategies)){

    if(strategies[s] == "NP"){
      variance_fraction <- 0
    }

    else if(strategies[s] == "P"){

      if(is.null(varfrac)){
        stop("Error: 0 < varfrac < 1")
      }
      variance_fraction <- as.numeric(varfrac)
    }

    else if(strategies[s] == "2-step"){
      variance_fraction <- 2
    }


    MixedImputationApproaches <- list("mice_phylo", list(data, nbrMI,
                                                         variance_fraction, tree, hint),
                                      "missForest_phylo", list(data, variance_fraction,
                                                               maxiter = 10, ntree = 100,
                                                               mtry = floor(ncol(data)/3), tree, hint),
                                      "kNN_phylo", list(data, k, numFun, catFun,
                                                        variance_fraction, tree, hint),
                                      "gain_phylo", list(data, variance_fraction, tree,
                                                    batch_size = round(ncol(data)*0.2),
                                                    hint_rate = 0.9, alpha = 100, epochs = 10000, hint))

    #to use only the imputation methods in ImputationApproachesNames
    methodsIndex <- which(MixedImputationApproaches %in% ImputationApproachesNames)

    for(method in methodsIndex){
      #univariate, missrandomForest and imputeKNN don't work.
      if((MixedImputationApproaches[[method]] == "mice_phylo" |
          MixedImputationApproaches[[method]] == "missForest_phylo" |
          MixedImputationApproaches[[method]] == "kNN_phylo") & ncol(data) == 1 & variance_fraction == 0){
        next
      }

      imputeName <- MixedImputationApproaches[[method]]
      imputeName[stringr::str_detect(imputeName, "mice")] <- "MICE"
      imputeName[stringr::str_detect(imputeName, "miss")] <- "MissForest"
      imputeName[stringr::str_detect(imputeName, "kNN")] <- "KNN"
      imputeName[stringr::str_detect(imputeName, "gain")] <- "GAIN"


      imputedValue <- do.call(MixedImputationApproaches[[method]], MixedImputationApproaches[[method + 1]])$imputedData
      # if(length(contiColumns) != 0){
      #   imputedValue[, contiColumns] <- exp(imputedValue[, contiColumns])
      # }
      names(imputedValue) <- colNames
      imputedData <- c(imputedData, list(imputedValue))

      #change names



      imputedName <- paste(imputeName, strategies[s], sep = "_")
      imputedNames <- c(imputedNames, imputedName)
    }
  }

  names(imputedData) <- imputedNames

  if(!is.null(save)){
    save(imputedData, file = paste0(save, ".RData"))
  }

  if(is.null(save)){
    return(imputedData)
  }

}