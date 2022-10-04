#' @title Imputation of missing data in continuous traits
#' @description This function imputes missing data for continuous traits applying the phylopars approach.
#' @usage pi_continuous_traits(missingData, tree)
#' @param missingData data.frame of 1 or more numeric columns containing NAs
#' @param tree phylo object
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the
#' imputation
#' @export
pi_continuous_traits <- function(missingData, tree){

  Rphylopars <- tryCatch(
    {

      if(length(setdiff(rownames(missingData), tree$tip.label)) != 0){
        rownames(missingData) <- tree$tip.label
      }

      #add species (tips) in the dataframe
      missingData <- tibble::add_column(missingData, species = row.names(missingData) , .before = 1)

      #impute
      models <- c("BM", "OU")
      AICs <- c()
      imputations <- list()

      phylo_correlated <- TRUE
      pheno_correlated <- FALSE

      for(i in 1:length(models)){

        imputeData <- Rphylopars::phylopars(trait_data = missingData, tree = tree, model = models[i],
                                            phylo_correlated = phylo_correlated, pheno_correlated = pheno_correlated)

        imputations[[i]] <- imputeData

        #Calculate AIC
        LogLik <- imputeData$logLik
        K <- imputeData$npars #get the number of degree of freedom (= the nbr of parameter in the model)
        AIC <- 2 * K - 2 * LogLik
        AICs <- c(AICs, AIC)
      }

      #keep only the imputed data
      imputedData <- imputations[[which.min(AICs)]]$anc_recon

      #keep only the tip labels and not the node labels
      imputedData <- imputedData[1:nrow(missingData), , drop = FALSE]

      #put again the tips names
      rownames(imputedData) <- tree$tip.label

      parameters <- list(model = imputations[[which.min(AICs)]]$model)

      list(imputedData = imputedData, parameters = parameters)
    },
    error = function(cond){
      message("Error running phylopars")
      list(imputedData = missingData, parameters = "Error in Cholesky decomposition")
    }
  )
  return(Rphylopars)
}
