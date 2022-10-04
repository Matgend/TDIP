#' @title Empirical data and phylogeny in R list
#'
#' @description This function saves an empirical dataset and a phylogeny in a list of the same structure like when the data
#' is simulated and converts the character and integer columns as factors.
#'
#' @usage data_preprocessing(empData, empTree = NULL)
#'
#' @param empData a table of class data.frame. Must have a column named Species
#' @param empTree a phylogenetic tree, by default empTree = NULL.
#' @param save path to save the data
#' @return A list which mimic the structure of the simulated data list. Contain at least the empirical data as data
#' .frame and the phylogenetic tree as phylo object. The continuous traits are log transformed.
#' @export
data_preprocessing <- function(empData, empTree = NULL, save = NULL){

  empData <- as.data.frame(empData)

  colNames <- names(empData)

  TreeList <- NULL
  if(!is.null(empTree)){

    #class the data.frame in the same order than the phylogenetic tips
    empData <- empData[match(empTree$tip.label, empData$Species),]

    #convert tree as ultrametric (needed for PI methods)
    if(!ape::is.ultrametric(empTree)){
      empTree <- phytools::force.ultrametric(empTree)
    }

    #reorder dataframe like this the order of the species match with the one the tips
    empData[match(empTree$tip.label, empData$Species), ]

    TreeList <- list('0' = empTree)
  }

  #convert character and integer variables as factors
  classVar <- lapply(empData, class)
  indexInFactor <- which(classVar == "integer" | classVar == "character" | classVar == "factor")
  empData[indexInFactor] <- lapply(empData[indexInFactor], factor)

  indexInFactor

  #log transform the continuous traits
  indexConti <- which(!1:ncol(empData) %in% indexInFactor)

  #check if negative values in log trait
  indexNega <- which(colnames(empData) %in% names(colSums(empData[indexConti] < 0) > 0))

  #remove trait with negative values for log transformation
  indexConti <- which(!indexConti %in% indexNega)
  empData[indexConti] <- log(empData[indexConti])

  row.names(empData) <- empData$Species

  Data <- list(FinalData = empData, TreeList = TreeList)

  if(!is.null(save)){
    save(Data, file = paste0(save, ".RData"))
  }
  return(Data)
}
