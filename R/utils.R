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
#' @usage corMatrixConti(discreteTrait, highCor)
#'
#' @param discreteTrait vector of integers
#' @param highCor numerical, strength of the correlation.
#' @return A vector of continuous variables normally distributed and correlated to the discrete trait.
corMatrixConti <- function(discreteTrait, highCor){

  discreteTrait <- as.numeric(discreteTrait)
  discreteTrait2 <- discreteTrait
  #shuffle discrete trait (randomize order of the states)
  nbrState <- length(unique(discreteTrait))
  S <- sample(0:(nbrState - 1), nbrState)
  for(i in 1:nbrState){
    discreteTrait[discreteTrait2 == S[i]] <- (i - 1)
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
#' @usage corContiTraitsOneTrait(nbrTraits, discreteTrait, highCor)
#'
#' @param nbrTraits vector, defining the number of traits (columns)
#' @param discreteTraits vector from which the simulated traits are correlated
#' @param highCor numerical, strength of the correlation.
#' @return A data frame of continuous traits where each column is a trait. Continuous traits are numeric.
corContiTraitsOneTrait <- function(nbrTraits, discreteTrait, highCor){

  #convert discreteTrait
  discreteTrait <- as.numeric(discreteTrait)

  contiMatrix <- lapply(1:nbrTraits, function(x) corMatrixConti(discreteTrait, highCor))
  contiMatrix <- matrix(unlist(contiMatrix), ncol = nbrTraits)
  #colnames(contiMatrix) <- paste0("corC", 1:ncol(contiMatrix))

  #get correlation matrix of final matrix
  finalMatrixConverted <- apply(contiMatrix, 2, as.numeric)
  correlationMatrix = cor(finalMatrixConverted)

  param <- list(correlationMatrix = correlationMatrix)

  return(list(data = as.data.frame(contiMatrix, stringsAsFactors=TRUE), parameters = param))
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
#' @return A list containing a data frame in which each discrete variable is converted as one hot encoding and vector of
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
#' @return A data.frame in which each variable that were one hot encoding are now factors.
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
