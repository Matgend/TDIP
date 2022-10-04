#'@title Impute NaNs in a simulated data according the PhyloNa mechanism
#'
#'@description Simulate NaNs in a complete data composed of discrete and continuous traits which are correlated
#'or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism PhyloNa. It used the function
#'getTipsNA().
#'
#'@usage phyloNa_miss_meca(missingRate, ds, tree)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data. Should be 0 < missing rate <= 0.5
#'@param ds dataframe in which NA should be created
#'@param tree object of class phylo
#'@return A dataset with NA following a pattern of MNAR.
#'@export
phyloNa_miss_meca <- function(missingRate, ds, tree){

  #in case missingRates larger than 50%
  if(missingRate > 0.5){
    stop("Can't simulate PhyloNaN")
  }

  #tips to include NaN
  tips <- getTipsNA(tree, missingRate)

  nameTrees <- paste("PhyloNaN", length(tips), round(missingRate, 2), sep = "/")

  ds[tips, ] <- NA

  PhyloNa <- list(ds)
  names(PhyloNa) <- nameTrees

  return(PhyloNa)
}
