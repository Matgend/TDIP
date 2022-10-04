#'@title Impute NaNs in a empirical data according to MCAR, MAR, MNAR and PhyloNa
#'
#'@description Simulate NaNs in a complete data set composed of continuous and discrete traits. The Nas are imputed according a
#'missing rate and respecting 4 categories of missing data, MCAR, MAR, MNAR and PhyloNa.
#'
#'@usage na_insertion(missingRates, dataset, missingTraits, MARtrait, save = NULL)
#'
#'@param missingRate float corresponding to the rate of missing value to introduce in the data
#'@param dataset data.frame
#'@param missingTraits numerical, total number of traits in which missing data will be generated
#'@param MARTraits index or vector, variables in which the missing values will be generated. (only for MAR) same length
#'than MARctrlTraits
#'@param MARctrlTraits index or vector, variables on which the missing values will be designed. (only for MAR) same length
#'than MARTraits
#'@param traitsNoNA integer mentioning the columns in which we don't want NA generation.
#'@param tree phylo object, by default NULL (in case no phylogenetic tree)
#'@param save character correspond to the name of the saved file in .RData format
#'@return a nested list composed of the partitioned data with the 4 categories of missing data (MCRA, MAR, MNAR and PhyloNa
#') according to a precise missing rate and of list of the parameters used for the missing values generation
#'@export
na_insertion <- function(missingRate, dataset, missingTraits, MARTraits, MARctrlTraits = NULL, traitsNoNA = NULL,
                             tree = NULL, save = NULL){

  colList <- colnames(dataset)

  #remove traits in which we don't want NA inside
  if(!is.null(traitsNoNA)){
    namesNoNA <- colList[traitsNoNA]
    colList <- colList[! colList %in% namesNoNA]
  }

  if (missingTraits > length(colList)){
    missTraits <- length(colList)
  }

  if(missingTraits <= length(colList)){
    missTraits <- missingTraits
  }

  if(length(colList) == missTraits){
    colMis <- colList
  }

  if(length(colList) < missTraits){
    colMis <- sample(colList, missTraits)
  }

  #MCAR
  missingMCAR <- mcar_miss_meca(missingRate, dataset, colMis)
  MCAR <- list(missingMCAR)
  nameMCAR <- paste("MCAR", length(colList),
                    round(missingRate, 2) ,sep = "/")
  names(MCAR) <- nameMCAR

  #MNAR
  missingMNAR <- mnar_miss_meca(missingRate, dataset, colMis)
  MNAR <- list(missingMNAR)
  nameMNAR <- paste("MNAR", length(colList),
                    round(missingRate, 2) ,sep = "/")
  names(MNAR) <- nameMNAR


  #MAR

  if(is.character(MARTraits)){
    cols_misIndex <- match(MARTraits, colList)
  }
  if(is.numeric(MARTraits)){
    cols_misIndex <-MARTraits
  }


  if(!is.null(MARctrlTraits)){

    if(length(MARTraits) != length(MARctrlTraits)){
      stop("MARTraits and MARctrlTraits must have the same size")
    }

    else if(is.character(MARctrlTraits)){
      cols_ctrlIndex <- match(MARctrlTraits, colList)
    }

    else if(is.numeric(MARctrlTraits)){
      cols_ctrlIndex <- MARctrlTraits
    }
  }

  if(is.null(MARctrlTraits)){

    cols_ctrlIndex <- c()

    for(c in 1:length(cols_misIndex)){
      pattern <- stringr::str_extract(names(dataset)[cols_misIndex[c]], "\\.\\d")
      MARctrlindex <- grep(pattern, names(dataset))
      MARctrlindex <- MARctrlindex[-cols_misIndex[c]]
      MARctrlindex <- MARctrlindex[which(!(MARctrlindex %in% cols_ctrlIndex ))]


      if(length(MARctrlindex) == 0){
        stop("MARctrlTraits and MARTraits must have the same length and the sum of both should
             be equal to number of trait being correlated together")
      }

      ctrl_index <- sample(MARctrlindex, length(cols_misIndex[c]))
      cols_ctrlIndex <- c(cols_ctrlIndex, ctrl_index)
    }

    cols_ctrlIndex <- sort(cols_ctrlIndex)
  }

  missingMAR <- mar_miss_meca(missingRate, dataset, cols_misIndex, cols_ctrlIndex)
  MAR <- list(missingMAR)
  nameMAR <- paste("MAR", length(colList),
                   round(missingRate, 2) ,sep = "/")
  names(MAR) <- nameMAR

  #PhyloNa
  if(!is.null(tree)){

    PhyloNaN <- phyloNa_miss_meca(missingRate, dataset, tree)

    #Data with NA imputed
    DataNaN <- list(MCAR = MCAR, MAR = MAR, MNAR = MNAR, PhyloNaN = PhyloNaN)
  }

  if(is.null(tree)){
    #Data with NA imputed
    DataNaN <- list(MCAR = MCAR, MAR = MAR, MNAR = MNAR)
  }

  NaNImputed <- list(DataNaN = DataNaN)

  #Save data
  ##########
  if(!is.null(save)){
    save(NaNImputed, file = paste0(save, ".RData"))
  }

  return(NaNImputed)
}
