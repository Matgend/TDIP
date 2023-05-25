######################################
# Function used in phyloNa_miss_meca()
######################################

#' @title Identify species for phylogenetic NAs
#'
#' @description This function returns the names of the species that should get NAs for their traits
#'
#' @usage getTipsNA(Tree, MinTips)
#'
#' @param Tree Phylogeny
#' @param missingRates numerical vector corresponding to the rate of missing value to introduce in the data
#' @return A vector containing the names of the species that should get NAs for their traits
#'
getTipsNA <- function (Tree, missingRates){
  Nodes <- ape::Nnode(Tree)
  Tips <- ape::Ntip(Tree)
  MaxNodeID <- Tips + Nodes
  NodeIDs <- (Tips + 1):MaxNodeID
  MinTips <- round(missingRates * Tips, 1)
  EnoughTips <- FALSE
  while (!EnoughTips) {
    FocalNode <- sample(NodeIDs, 1)
    Desc <- phangorn::Descendants(Tree, FocalNode, type = "tips")[[1]]
    SampledTips <- Tree$tip.label[Desc]
    if(length(SampledTips) >= MinTips){
      EnoughTips <- TRUE
    }
  }
  if(length(SampledTips) > MinTips){
    SampledTips <- sample(SampledTips, MinTips)
  }
  return(SampledTips)
}


#####################################
# Functions for data simulation
#####################################
#' @title Simulate variance-covariance matrix for morphological evolution
#'
#' @description This function generates a variance-covariance matrix for
#' morphological evolution, defining the rate at which traits evolve and their
#' evolutionary correlation.
#'
#' @usage simSigma(Ntraits, Cor = NULL, Sigma2 = NULL, uncovTraits = NULL, FracNocov = NULL)
#'
#' @param Ntraits number of traits
#' @param Cor correlation between traits. Default between -1 and 1.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' @param Sigma2 Brownian motion rate. Default between 1e-4 and 0.5.
#' Optional, can be fixed to be equal for all traits or one value per trait
#' @param uncovTraits number of traits having a covariance equal to 0 with the others.
#' @param FracNocov fraction of covariance being zero
#'
#' @return A matrix Ntrait * Ntrait for simulating trait evolution
simSigma <- function(Ntraits, Cor = NULL, Sigma2 = NULL, uncovTraits = NULL, FracNocov = NULL){
  if(!is.null(Sigma2)){
    if(length(Sigma2) != Ntraits && length(Sigma2) != 1){
      stop("Sigma2 should be of length 1 or Ntraits")
    }else if(length(Sigma2) == 1) {
      Sigma2 <- rep(Sigma2, Ntraits)
    }
  }
  else{
    Sigma2 <- runif(Ntraits, min = 1e-4, max = 0.5)
  }

  Sigmas <- matrix(Sigma2, nrow = 1)
  if(Ntraits > 1) {
    Cov <- matrix(1, ncol = Ntraits, nrow = Ntraits)
    Q <- Ntraits*(Ntraits-1)/2 #for me it's +1 and not -1
    if(!is.null(Cor)) {
      if(length(Cor) != Q && length(Cor) != 1) {
        stop("Correlation among traits should be of length 1 or Ntraits*(Ntraits-1)/2")
      }
      SimCov <- Cor
    }
    else{
      SimCov <- runif(Q, min = -1, max = 1) # Trait correlation
    }
    Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
    Cov[upper.tri(Cov, diag = FALSE)] <- SimCov
    Sigmas <- diag(Sigma2)  %*% Cov  %*% diag(Sigma2) # Correlation to covariance

    # Force variance-covariance to be positive definite # can I have some precision?
    Tol <- 1e-6
    Ev <- eigen(Sigmas, symmetric = TRUE)$values
    if(!all( Ev >= -Tol * abs(Ev[1L]))) {
      Sigmas <- as.matrix(Matrix::nearPD(Sigmas)$mat)
    }
  }

  if(!is.null(uncovTraits) && uncovTraits < 0){
    stop("UncovTraits should be equal or greater than 0")
  }

  if(!is.null(uncovTraits) && uncovTraits %% 1 != 0){
    stop("UncovTraits should be an integer")
  }

  if(length(uncovTraits) == 1 && uncovTraits == 0){
    uncovTraits = NULL
  }

  if(!is.null(uncovTraits) && uncovTraits > Ntraits){
    stop("Uncovtraits should be equal or smaller than Ntraits")
  }

  if(!is.null(uncovTraits) && Ntraits > 1) {
    columns <- 1:uncovTraits
    for(i in columns){
      valueToKeep <- Sigmas[i, i]
      Sigmas[i,] <- 0
      Sigmas[,i] <- 0
      Sigmas[i, i] <- valueToKeep
    }
  }

  if(length(FracNocov) == 1 && FracNocov == 0){
    FracNocov = NULL
  }

  if(!is.null(FracNocov) && Ntraits > 1) {
    Q <- Ntraits*(Ntraits-1)/2
    Nzero <- round(Q * FracNocov) #to have the number of variable having no covariance.
    Mask <- matrix(1, nrow = Ntraits, ncol = Ntraits)
    Lt <- which(lower.tri(Mask))
    SetZero <- sample(Lt, Nzero)
    Mask[SetZero] <- 0
    Mask <- t(Mask)
    Mask[SetZero] <- 0
    Sigmas <- Sigmas * Mask
  }
  return(Sigmas)
}

#' @title Simulate alpha matrix for morphological evolution
#'
#' @description This function generates an alpha matrix for
#' morphological evolution, defining the strength of attraction
#' towards a morphological optimum
#'
#' @usage simSigma(Ntraits, alpha = NULL)
#'
#' @param Ntraits number of traits
#' @param alpha strength of attraction towards morphological optimum. Default between 0.5 and 2.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#'
#' @return A matrix Ntrait * Ntrait for simulating trait evolution
#'
simAlpha <- function(Ntraits, alpha = NULL) {
  if(is.null(alpha)){
    alpha <- exp(runif(Ntraits, log(0.5), log(2)))
    TreeHeight <- 1 # In case we need to change the age of the phylogenies
    alpha <- alpha * log(2) * 1/TreeHeight
    AlphaMat <- diag(alpha, nrow = Ntraits, ncol = Ntraits)
  }
  else {
    if(length(alpha) != Ntraits && length(alpha) != 1) {
      stop("alpha should be of length 1 or Ntraits")
    }
    if(any(alpha < 0)) {
      stop("Alpha should be larger or equal to 0")
    }
    if(length(alpha) == 1) {
      AlphaMat <- diag(rep(alpha, Ntraits), nrow = Ntraits, ncol = Ntraits)
    }
    else{
      AlphaMat <- diag(alpha, nrow = Ntraits, ncol = Ntraits)
    }
  }
  return(AlphaMat)
}

#' @title Simulate discrete traits for all species according to the MK model
#'
#' @description This function generates matrix of discrete values (states) of
#' morphological evolution. Define also a matrix of nodes
#'
#' @usage simDiscreteTraits(Ntraits, Nstates, rate_model, max_rate, SimTree)
#'
#' @param Ntraits number of traits
#' @param Nstates number of states of the defined traits (could be a vector or a value)
#' @param rate_model vector, rate model that the transition matrix must satisfy. Chose between "ER" (=all permitted
#' transitions occur at the same rate), "SYM" (=hat backward & forward transitions occur at the same rate) and "ARD" (=all
#' allowed transitions can occur at different rates).
#' @param tree stochastic birth-death trees
#' @param equal if TRUE, the transition matrix is equal for all the traits. If FALSE all the traits have a different
#' transition matrix. By default is TRUE.
#' @param Ordinal simulate ordinal data. Default is FALSE
#' @return A list containing morphological evolution for each trait and each species, a matrix Ntrait x Ntrait for simulating
#' trait evolution and sigma matrix Ntraits x Ntraits
simDiscreteTraits <- function(Ntraits, Nstates, rate_model, max_rate, tree, equal = TRUE, Ordinal = FALSE){

  if(rate_model != "ER" & rate_model != "SYM" & rate_model != "ARD"){
    stop("The rate model should be one of the following: ER, SYM or ARD")
  }

  if(length(rate_model) != 1 && length(rate_model) != Ntraits){
    stop("rate_model should be of length 1 or Nstates")
  }

  if(length(rate_model) == 1 && length(rate_model) != Ntraits) {
    rate_model <- rep(rate_model, Ntraits)
  }

  if(length(Nstates) != 1 && length(Nstates) != Ntraits){
    stop("Nstates should be of length 1 or Nstates")
  }

  if(length(Nstates) == 1 && length(Nstates) != Ntraits) {
    Nstates <- rep(Nstates, Ntraits)
  }

  tip_mat <- matrix(0, nrow = length(tree$tip.label), ncol = Ntraits)
  node_mat <- matrix(0, nrow = tree$Nnode, ncol = Ntraits)
  for(i in 1:Ntraits) {
    # define transition matrix
    if(equal){
      set.seed(1)
      Q <- castor::get_random_mk_transition_matrix(Nstates[i], rate_model = rate_model[i], max_rate = max_rate)
    }
    else{
      Q <- castor::get_random_mk_transition_matrix(Nstates[i], rate_model = rate_model[i], max_rate = max_rate)
    }
    if(Ordinal){ #going to be have more rate of change which are equal to 0.
      for (y in 2:Nstates[i]) {
        Q[row(Q) == (col(Q) - y)] <- 0
        Q[row(Q) == (col(Q) + y)] <- 0
      }
    }
    tip_states <- castor::simulate_mk_model(tree, Q)
    tip_mat[, i] <- tip_states$tip_states
    node_mat[, i] <- tip_states$node_states # Do we need the ancestral states?
  }

  #get state starting to 0
  tip_mat <- tip_mat - 1
  node_mat <- node_mat - 1

  return(list(tip_mat = tip_mat, node_mat = node_mat))
}

#' @title Converts a single continuous variable in discrete variable
#'
#' @description This function converts a continuous variable in discrete variable. The variable can be converted as interval,
#' ordinal or nominal.
#'
#' @usage ConvertContinousInDiscreteValues(values, Nstates, subclass)
#'
#' @param values vector of float
#' @param Nstates number of states of the defined traits
#' @param subclass intervals: states are fairly split. ordinal: ordered states,
#' split is random. non_eq_nominal: split is random, no order (shuffled)
#' @return A vector with integer instead of floats
ConvertContinousInDiscreteValues <- function(values, Nstates, subclass){

  if(subclass != "non_eq_nominal" & subclass != "ordinal" & subclass != "interval" & subclass != "eq_nominal"){
    stop("The subclass should be one of the following: nominal, ordinal or interval")
  }

  if(subclass == "interval"){
    breaks <- seq(min(values), max(values), length.out = Nstates + 1)
    conversion <- as.character(findInterval(values, breaks[-c(1, length(breaks))]))
  }

  if(subclass == "ordinal"){
    breaks <- runif((Nstates - 1), min(values), max(values))
    conversion <- as.character(findInterval(values, sort(breaks)))
  }

  if(subclass == "non_eq_nominal"){
    breaks <- sample(sort(values), Nstates-1)
    nominal_values <- findInterval(values, sort(breaks))
    shuffle <- sample(0:(Nstates-1), Nstates, replace = FALSE)
    conversion <- as.character(1:length(nominal_values)) # vector with shuffling
    for (i in 1:length(unique(nominal_values))){
      conversion[nominal_values == unique(nominal_values)[i]] <- shuffle[i]
    }
  }

  return(conversion)
}

#' @title Convert several continuous variables in discrete variables
#'
#' @description This function converts continuous variables in discrete variables. The variables can be converted as
#' interval, ordinal or nominal.
#'
#' @usage ChangeContinuousTraitInDiscrete(Matrix, columnsIndex, Nstates, subclass)
#'
#' @param Matrix array matrix of continuous variables
#' @param columnsIndex vector of integers mentioning the columns to convert in integers
#' @param Nstates number of states of the defined traits (could be a vector or a value)
#' @param subclass intervals: states are fairly split. ordinal: ordered states,
#' split is random. non_eq_nominal: split is random, no order (shuffled)
#' @return A data frame with the integer variables and continuous variables.
ChangeContinuousTraitInDiscrete <- function(Matrix, columnsIndex, Nstates, subclass){

  if(length(Nstates) >= 1 && length(Nstates) < length(columnsIndex)){
    stop("Nstates length should be equal to 1 or equal to the number of columnsIndex")
  }

  if(length(Nstates) == 1){
    Nstates <- rep(Nstates, length(columnsIndex))
  }

  if(length(subclass) == 1){
    subclass <- rep(subclass, length(columnsIndex))
  }

  #convert matrix in dataframe
  dataframe <- as.data.frame(Matrix)

  for(ci in 1:length(columnsIndex)) {
    dataframe[, columnsIndex[ci]] <- ConvertContinousInDiscreteValues(
      dataframe[, columnsIndex[ci]], Nstates[ci], subclass[ci])
  }
  return(dataframe)
}

#' @title Probability matrix
#'
#' @description This function generates a random array of size nbrState * nbrState for which the sum of each row is equal
#' to 1.
#'
#' @usage corMatrixDisc(nbrState, highCor)
#'
#' @param nbrStates integer, providing the number of states to include in the simulated vectors.
#' @param highCor numerical, strength of the correlation. For instance, if nbrState = 3 and highCor = 0.8, the other 2 traits will have
#' a correlation of (1-0.8)/(length(states) - 1)
#' @param states vector of the states (unique) present in the discrete vector
#' @return An array of size nbrState * length(states) for which the sum of each row is equal to 1.
corMatrixDisc <- function(nbrState, highCor){
  lowCor <- (1 - highCor) / (nbrState - 1)
  ProbMat <- matrix(rep(lowCor, nbrState * nbrState), nbrState, nbrState)
  diag(ProbMat) <- highCor
  IdxReorder <- sample(1:nbrState, nbrState)
  ProbMat <- ProbMat[IdxReorder, ]
  colnames(ProbMat) <- 0:(nbrState - 1)
  rownames(ProbMat) <- 0:(nbrState - 1)
  return(ProbMat)
}

#' @title Correlated continuous variable
#'
#' @description This function generates a continuous trait which is correlated to a discrete trait according to a
#' correlation factor. The function applies the function rnorm_pre from the R package faux.
#'
#' @usage corMatrixConti(discreteTrait, highCor, shuffle)
#'
#' @param discreteTrait vector of integers
#' @param highCor numerical, strength of the correlation.
#' @param shuffle boolean, if true discrete states are shuffled, False, discrete states are not shuffled
#' @return A vector of continuous variables normally distributed and correlated to the discrete trait.
corMatrixConti <- function(discreteTrait, highCor, shuffle = TRUE){

  discreteTrait <- as.numeric(discreteTrait)
  discreteTrait2 <- discreteTrait
  #shuffle discrete trait (randomize order of the states)
  nbrState <- length(unique(discreteTrait))
  if(shuffle){
    S <- sample(0:(nbrState - 1), nbrState)
    for(i in 1:nbrState){
      discreteTrait[discreteTrait2 == S[i]] <- (i - 1)
    }
  }
  contiTrait <- faux::rnorm_pre(discreteTrait,
                                mu =  mean(discreteTrait),
                                sd = sd(discreteTrait),
                                r = highCor,
                                empirical = FALSE,
                                threshold = 1e-12)

  return(contiTrait)
}


#' @title Discrete traits correlated to a single trait discrete trait
#'
#' @description This function generates a data frame of discrete traits
#'  which are all correlated with a single traits but uncorrelated with each other.
#'
#' @usage corDiscTraitsOneTrait(nbrTraits, discreteTrait, highCor)
#'
#' @param nbrTraits vector, defining the number of traits (columns)
#' @param discreteTraits vector from which the simulated traits are correlated
#' @param highCor numerical, strength of the correlation. For instance, if nbrState = 3 and highCor = 0.8, the other 2 traits will have
#'  a correlation of (1-0.8)/(nbrTraits - 1)
#' @return A data frame of discrete traits where each column is a trait. Discrete traits are factors.
corDiscTraitsOneTrait <- function(nbrTraits, discreteTrait, highCor){

  #convert discreteTrait
  discreteTrait <- as.numeric(discreteTrait)

  nbrState <- length(unique(discreteTrait))
  ProbMat <- lapply(1:nbrTraits, function(x) corMatrixDisc(nbrState, highCor))
  l <- lapply(1:nbrTraits, function(x){
    sapply(1:length(discreteTrait), function(y)
      as.character(sample(0:(nbrState - 1), 1, prob = ProbMat[[x]][(discreteTrait[y] + 1), ])))
  })

  discMatrix <- matrix(unlist(l), ncol = nbrTraits)
  #colnames(finalMatrix) <- paste0("corD", 1:ncol(finalMatrix))
  param <- list(ProbMat)

  #get correlation matrix of final matrix
  finalMatrixConverted <- apply(discMatrix, 2, as.numeric)
  correlationMatrix = cor(finalMatrixConverted)

  param <- c(param, list(correlationMatrix = correlationMatrix))

  return(list(data = as.data.frame(discMatrix, stringsAsFactors=TRUE), parameters = param))
}


#' @title Continuous traits correlated to a single trait discrete trait
#'
#' @description This function generates a data frame of continuous traits
#'  which are all correlated with a single traits but uncorrelated with each other.
#'
#' @usage corContiTraitsOneTrait(nbrTraits, discreteTrait, highCor, shuffle)
#'
#' @param nbrTraits vector, defining the number of traits (columns)
#' @param discreteTraits vector from which the simulated traits are correlated
#' @param highCor numerical, strength of the correlation.
#' @param shuffle boolean, if true discrete states are shuffled, False, discrete states are not shuffled
#' @return A data frame of continuous traits where each column is a trait. Continuous traits are numeric.
corContiTraitsOneTrait <- function(nbrTraits, discreteTrait, highCor, shuffle = TRUE){

  #convert discreteTrait
  discreteTrait <- as.numeric(discreteTrait)

  contiMatrix <- lapply(1:nbrTraits, function(x) corMatrixConti(discreteTrait, highCor, shuffle))
  contiMatrix <- matrix(unlist(contiMatrix), ncol = nbrTraits)
  #colnames(contiMatrix) <- paste0("corC", 1:ncol(contiMatrix))

  #get correlation matrix of final matrix
  finalMatrixConverted <- apply(contiMatrix, 2, as.numeric)
  correlationMatrix = cor(finalMatrixConverted)

  param <- list(correlationMatrix = correlationMatrix)

  return(list(data = as.data.frame(contiMatrix, stringsAsFactors=TRUE), parameters = param))
}

#' @title Rescale a phylogenetic tree
#'
#' @description This function rescales a phylogenetic tree of class phylo according to kappa or lambda (Pagel)
#'
#' @usage rescaleTree(tree, subdata)
#'
#' @param tree phylo object
#' @param subdata a dataframe composed of the columns of the parameter called dataframe of the data_simulator function
#' @return a phylo object rescaled according to the model define in the data frame.
rescaleTree <- function(tree, subdata){

  #rescale phylogeny
  lambdaCheck <- mean(subdata$lambda)
  kappaCheck <- mean(subdata$kappa)

  if(lambdaCheck != 0 & kappaCheck != 0 & lambdaCheck != 1 & kappaCheck != 1){
    stop("lambda or kappa should be equal to 1")
  }

  else if(lambdaCheck != 1){
    subdataTree <- geiger::rescale(tree, "lambda", lambdaCheck)
  }

  else if (kappaCheck != 1){
    subdataTree <- geiger::rescale(tree, "kappa", kappaCheck)
  }

  else{
    subdataTree <- tree
  }

  return(subdataTree)
}


###################
#Functions for GAIN
###################

#' @title Conversion factors in one hot encoding
#'
#' @description This function converts discrete variables in one hot encoding using the function one_hot from the R package
#' mltools
#'
#' @usage generateOneHotEncoVariables(NaNData)
#'
#' @param NaNData data frame of one or several factors columns
#' @return A list containing a dataframe in which each discrete variable is converted as one hot encoding and vector of
#' characters.
generateOneHotEncoVariables <- function(NaNData){

  nativeColNames <- colnames(NaNData)
  nativeRowNames <- rownames(NaNData)

  #in case NaNData is a vector
  if(is.null(dim(NaNData))){
    stop("The column(s) shoud be of class data.frame")
  }

  #convert automatically the factor or character columns into one hot encoding
  oneHotdata <- mltools::one_hot(data.table::as.data.table(NaNData))
  oneHotdata <- as.data.frame(oneHotdata)
  row.names(oneHotdata) <- nativeRowNames

  return(list(data = oneHotdata, nativeColNames = nativeColNames))
}

#' @title Conversion one hot encoding in factors
#'
#' @description This function converts data frame of one hot encoding
#' representing a categorical variable.
#'
#' @usage convertOneHotEncoIntoFactors(oneHotData, nativeColNames)
#'
#' @param oneHotData data frame of one one hot encoding
#' @param nativeColNames columns names of the original column names of the one hot encoding(before the conversion in one hot
#' encoding).
#' @return A dataframe in which each variable that were one hot encoding are now factors.
convertOneHotEncoIntoFactors <- function(oneHotData, nativeColNames){

  #isolate one hot encoding variables
  colNames <- colnames(oneHotData)

  oneHotColNames <- colNames[grep("\\_", colNames)]

  #select the corresponding columns in the nativeColNames vector
  factorColumn <- nativeColNames[nativeColNames %in% unique(stringr::str_extract(oneHotColNames, "^.*(?=\\_)"))]

  factorDataframe <- data.frame(matrix(NA, nrow = nrow(oneHotData),
                                       ncol = length(nativeColNames)))
  names(factorDataframe) <- nativeColNames

  for(c in 1:length(nativeColNames)){

    if(nativeColNames[c] %in% unique(stringr::str_extract(oneHotColNames, "^.*(?=\\_)"))){

      #select data by group of one hot encoding
      subgroupOneHotColNames <- oneHotColNames[grep(nativeColNames[c], oneHotColNames)]

      #isolate states of the trait
      states <- stringr::str_extract(subgroupOneHotColNames, "\\d+$")

      discreteTrait <- c()
      #go through the rows
      for(r in 1:nrow(oneHotData)){

        valueState <- states[which(unlist(oneHotData[r, subgroupOneHotColNames]) == 1)]

        #in case 2 or 3 one in the row
        if(length(valueState) > 1){
          valueState <- sample(valueState, 1)
        }

        #in case only 0 in the row
        else if(length(valueState) == 0){
          valueState <- sample(states, 1)
        }

        discreteTrait <- c(discreteTrait, valueState)
      }

      factorDataframe[ , nativeColNames[c]] <- as.factor(discreteTrait)

    }
    else{
      factorDataframe[ , nativeColNames[c]] <- oneHotData[, nativeColNames[c]]
    }
  }

  return(factorDataframe)

}

#modified missForest function from missForest package written by D.Stekhoven
##
## MissForest - nonparametric missing value imputation for mixed-type data
##
## This R script contains the actual missForest function.
##
## Author: D.Stekhoven, stekhoven@nexus.ethz.ch
##
## Acknowledgement: Steve Weston for input regarding parallel execution (2012)
##############################################################################

myMissForest <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                       xtrue = NA, parallelize = c('no', 'variables', 'forests'))
{ ## ----------------------------------------------------------------------
  ## Arguments:
  ## xmis         = data matrix with missing values
  ## maxiter      = stop after how many iterations (default = 10)
  ## ntree        = how many trees are grown in the forest (default = 100)
  ## variablewise = (boolean) return OOB errors for each variable separately
  ## decreasing   = (boolean) if TRUE the columns are sorted with decreasing
  ##                amount of missing values
  ## verbose      = (boolean) if TRUE then missForest returns error estimates,
  ##                runtime and if available true error during iterations
  ## mtry         = how many variables should be tried randomly at each node
  ## replace      = (boolean) if TRUE bootstrap sampling (with replacements)
  ##                is performed, else subsampling (without replacements)
  ## classwt      = list of priors of the classes in the categorical variables
  ## cutoff       = list of class cutoffs for each categorical variable
  ## strata       = list of (factor) variables used for stratified sampling
  ## sampsize     = list of size(s) of sample to draw
  ## nodesize     = minimum size of terminal nodes, vector of length 2, with
  ##                number for continuous variables in the first entry and
  ##                number for categorical variables in the second entry
  ## maxnodes     = maximum number of terminal nodes for individual trees
  ## xtrue        = complete data matrix
  ##
  ## ----------------------------------------------------------------------
  ## Author: Daniel Stekhoven, stekhoven@nexus.ethz.ch

  ## stop in case of wrong inputs passed to randomForest
  n <- nrow(xmis)
  p <- ncol(xmis)
  if (!is.null(classwt))
    stopifnot(length(classwt) == p, typeof(classwt) == 'list')
  if (!is.null(cutoff))
    stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
  if (!is.null(strata))
    stopifnot(length(strata) == p, typeof(strata) == 'list')
  if (!is.null(nodesize))
    stopifnot(length(nodesize) == 2)

  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  }

  ## return feedback on parallelization setup
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c('variables', 'forests')) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    } else if (verbose) {
      if (parallelize == 'variables') {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      } else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p){
      stop("The number of parallel cores should not exceed the number of variables (p=", p, ")")
    }
  }

  ## perform initial S.W.A.G. on xmis (mean imputation)
  ximp <- xmis
  varType <- character(p)
  for (t.co in 1:p) {
    if (is.numeric(xmis[[t.co]])) {
      varType[t.co] <- 'numeric'
      ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
      next()
    }
    if (is.factor(xmis[[t.co]])) {
      varType[t.co] <- 'factor'
      ## take the level which is more 'likely' (majority vote)
      max.level <- max(table(ximp[[t.co]]))
      ## if there are several classes which are major, sample one at random
      class.assign <- sample(names(which(max.level == summary(ximp[[t.co]]))), 1)
      ## it shouldn't be the NA class
      if (class.assign != "NA's") {
        ximp[is.na(xmis[[t.co]]),t.co] <- class.assign
      } else {
        while (class.assign == "NA's") {
          class.assign <- sample(names(which(max.level ==
                                               summary(ximp[[t.co]]))), 1)
        }
        ximp[is.na(xmis[[t.co]]),t.co] <- class.assign
      }
      next()
    }
    stop(sprintf('column %s must be factor or numeric, is %s', names(xmis)[t.co], class(xmis[[t.co]])))
  }

  ## extract missingness pattern
  NAloc <- is.na(xmis)            # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
  if (decreasing)
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]

  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == 'variables') {
    '%cols%' <- get('%dorng%')
    idxList <- as.list(isplitVector(nzsort.j, chunkSize = getDoParWorkers()))
  }
  #   else {
  #     ## force column loop to be sequential
  #     '%cols%' <- get('%do%')
  #     idxList <- nzsort.j
  #   }

  ## output
  Ximp <- vector('list', maxiter)

  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType

  ## setup convergence variables w.r.t. variable types
  if (k == 1){
    if (unique(varType) == 'numeric'){
      names(convNew) <- c('numeric')
    } else {
      names(convNew) <- c('factor')
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c('numeric', 'factor')
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }

  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }

  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    if (iter != 0){
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    if (verbose){
      cat("  missForest iteration", iter+1, "in progress...")
    }
    t.start <- proc.time()
    ximp.old <- ximp

    if (parallelize == "variables"){
      for (idx in idxList) {
        results <- foreach(varInd = idx, .packages = 'randomForest') %cols% {
          obsi <- !NAloc[, varInd] # which i's are observed
          misi <- NAloc[, varInd] # which i's are missing
          obsY <- ximp[obsi, varInd] # training response
          obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
          misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
          typeY <- varType[varInd]
          if (typeY == 'numeric'){
            RF <- randomForest(
              x = obsX,
              y = obsY,
              ntree = ntree,
              mtry = mtry,
              replace = replace,
              sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
              nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
            ## record out-of-bag error
            oerr <- RF$mse[ntree]
            #           }
            ## predict missing values in column varInd
            misY <- predict(RF, misX)
          } else { # if Y is categorical
            obsY <- factor(obsY) ## remove empty classes
            summarY <- summary(obsY)
            if (length(summarY) == 1){ ## if there is only one level left
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            } else {
              RF <- randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                classwt = if (!is.null(classwt)) classwt[[varInd]] else
                  rep(1, nlevels(obsY)),
                cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                  rep(1/nlevels(obsY), nlevels(obsY)),
                strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                  if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              oerr <- RF$err.rate[[ntree,1]]
              #             }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)
            }
          }
          list(varInd = varInd, misY = misY, oerr = oerr)
        }
        ## update the master copy of the data
        for (res in results) {
          misi <- NAloc[,res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    } else { # if parallelize != "variables"
      for (s in 1 : p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd], drop = FALSE]
          misX <- ximp[misi, seq(1, p)[-varInd], drop = FALSE]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == 'forests') {
              xntree <- NULL
              RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()),
                            .combine = 'combine', .multicombine = TRUE,
                            .packages = 'randomForest') %dorng% {
                              randomForest( x = obsX,
                                            y = obsY,
                                            ntree = xntree,
                                            mtry = mtry,
                                            replace = replace,
                                            sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                              if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                            nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                            maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                            }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm = TRUE)
              #               OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- randomForest( x = obsX,
                                  y = obsY,
                                  ntree = ntree,
                                  mtry = mtry,
                                  replace = replace,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- predict(RF, misX)
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize == 'forests') {
                RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()),
                              .combine = 'combine', .multicombine = TRUE,
                              .packages = 'randomForest') %dorng% {
                                randomForest(
                                  x = obsX,
                                  y = obsY,
                                  ntree = xntree,
                                  mtry = mtry,
                                  replace = replace,
                                  classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                    rep(1, nlevels(obsY)),
                                  cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                    rep(1/nlevels(obsY), nlevels(obsY)),
                                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                              }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[! is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- randomForest(x = obsX,
                                   y = obsY,
                                   ntree = ntree,
                                   mtry = mtry,
                                   replace = replace,
                                   classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                     rep(1, nlevels(obsY)),
                                   cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                     rep(1 / nlevels(obsY), nlevels(obsY)),
                                   strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                   sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                     if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                   nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                   maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                ## record out-of-bag error
                OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
              }
              ## predict missing parts of Y
              misY <- predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    if (verbose){
      cat('done!\n')
    }

    iter <- iter + 1
    Ximp[[iter]] <- ximp

    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)){
      t.ind <- which(varType == t.type)
      if (t.type == 'numeric'){
        convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, t.ind])^2) / sum(ximp[, t.ind]^2)
      } else {
        dist <- sum(as.character(as.matrix(ximp[, t.ind])) != as.character(as.matrix(ximp.old[, t.ind])))
        convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
      }
      t.co2 <- t.co2 + 1
    }

    ## compute estimated imputation error
    if (!variablewise){
      NRMSE <- sqrt(mean(OOBerror[varType == 'numeric'])/
                      var(as.vector(as.matrix(xmis[, varType == 'numeric'])),
                          na.rm = TRUE))
      PFC <- mean(OOBerror[varType == 'factor'])
      if (k == 1){
        if (unique(varType) == 'numeric'){
          OOBerr <- NRMSE
          names(OOBerr) <- 'NRMSE'
        } else {
          OOBerr <- PFC
          names(OOBerr) <- 'PFC'
        }
      } else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c('NRMSE', 'PFC')
      }
    } else {
      OOBerr <- OOBerror
      names(OOBerr)[varType == 'numeric'] <- 'MSE'
      names(OOBerr)[varType == 'factor'] <- 'PFC'
    }

    if (any(!is.na(xtrue))){
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }

    ## return status output, if desired
    if (verbose){
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))){
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }#end while((convNew<convOld)&(iter<maxiter)){

  ## produce output w.r.t. stopping rule
  if (iter == maxiter){
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    } else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
    }
  } else {
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld,
                  error = suppressWarnings(mixError(Ximp[[iter - 1]], xmis, xtrue)))
    }
  }
  class(out) <- 'missForest'
  return(out)
}


