#' @title Imputation of missing data in discrete traits, mandatory phylogenetic information
#'
#' @description This function applies the function imputeOneDiscretTrait on several columns at one time.
#'
#' @usage pi_categorical_traits(missingData, tree)
#'
#' @param missingData data.frame of 1 or more factor columns containing NAs
#' @param tree phylo object
#' @return A data.frame of 1 or more factor columns with the NAs replaced by values, parameters used for the
#' imputation, matrix of probabilities of each states
#' @export
pi_categorical_traits <- function(missingData, tree){

  #select the columns with missing values
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))

  proba <- matrix(NA, nrow = nrow(missingData), ncol = ncol(missingData))
  parameters <- list()
  for(i in 1:length(names(NaNColIndex))){
    imputation <- imputeOneDiscreteTrait(missingData[, names(NaNColIndex)[i], drop = FALSE], tree)
    #print("ok")
    missingData[ ,names(NaNColIndex)[i]] <- imputation$imputedData

    if(i == 1){
      proba <- imputation$probabilities
    }

    else{
      proba <- cbind(proba, imputation$probabilities)
    }
    parameters <- c(parameters, imputation$parameters)
  }
  proba <- as.data.frame(proba, row.names = rownames(missingData))

  return(list(imputedData = missingData, parameters = parameters, probabilities = proba))
}


#' @title Imputation of missing data for one discrete trait
#'
#' @description This function imputes missing values of one vector according to a Markov chain model using the R function
#' corHMM from the corHMM package.
#'
#' @usage imputeOneDiscreteTrait(trait, Data)
#'
#' @param missingData data frame of 1 factor column containing NAs
#' @param tree phylo object
#' @return a data frame of 1 factor column with the NAs replaced by values.
#'
imputeOneDiscreteTrait <- function(missingData, tree){

  levelsTrait <- levels(missingData[,1])

  print(names(missingData))

  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(tree$tip.label, row.names(missingData))){

    #change order of the rows, match the order of the phylogeny
    missingData <- missingData[match(tree$tip.label, row.names(missingData)), drop = FALSE]
  }

  colName <- names(missingData)
  #if only one state represented in the trait
  if(sum(!is.na(unique(missingData[, 1]))) == 1){
    #print("one state")
    state <- missingData[which(!is.na(missingData)), ][1]
    missingData[which(is.na(missingData)), ] <- state
    return(list(imputedData = missingData, parameters = list("noModel")))
  }

  #add the tip names in the data frame
  missingData <- cbind(species = row.names(missingData), missingData)

  #convert missingData as character
  missingData[,2] <- as.character(missingData[,2])
  missingData[,2][which(is.na(missingData[,2]))] <- "?" #replace NA by "?"because corHMM don't like it
  #Define the rate model
  model <- "ER"
  FitCorHMM <- corHMM::corHMM(phy = tree, data = missingData, model = model, rate.mat = ,
                              rate.cat = 1, get.tip.states = TRUE)


  #Calculate AIC
  AIC <- FitCorHMM$AIC
  models <- c("SYM", "ARD")
  model <- "ER"
  for (i in 1:length(models)){

    FitCorHMMDiffModel <- corHMM::corHMM(phy = tree, data = missingData, model = models[i],
                                         rate.cat = 1, get.tip.states = TRUE)

    #Calculate AIC
    AICDiffModel <- FitCorHMMDiffModel$AIC
    if(AIC > AICDiffModel){
      AIC <- AICDiffModel
      FitCorHMM <- FitCorHMMDiffModel
      model <- models[i]
    }
  }

  #Imputation
  MostLikelyState <- apply(FitCorHMM$tip.states, 1, which.max)

  if(levelsTrait[1] == "0"){
    MostLikelyState <- as.data.frame(as.factor(MostLikelyState - 1)) #-1 because when converted a numeric the value are +1
  }

  if(levelsTrait[1] != "0"){
    MostLikelyState <- as.data.frame(as.factor(MostLikelyState)) #-1 because when converted a numeric the value are +1
  }

  colnames(MostLikelyState) <- colName

  #Parameters
  parameters <- list(model = model, rate.cat = 1)

  #print(MostLikelyState)
  return(list(imputedData = MostLikelyState, parameters = parameters, probabilities = FitCorHMM$tip.states))
}
