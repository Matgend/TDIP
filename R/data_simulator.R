#' @title Simulation of trait data
#'
#' @description This function simulates a trait data frame according several parameters.
#'
#' @usage data_simulator(param_tree, dataframe, save = NULL)
#'
#' @description The parameter `dataframe` is composed of 11 columns:
#' * nbr_traits: Number of traits simulated with specific parameters
#' * class: Type of traits, (continuous, non_eq_nominal(categorical) and ordinal)
#' * model: Evolutionary model (BM1, OU1, ARD, SYM, ER, Na)
#' * states: Number of states for discrete traits, (if continuous set it to 1)
#' * correlation: Index corresponding to the group of simulated traits which are correlated or not to other traits
#' * uncorr_traits: Among the "nbr_traits", it's the number of uncorrelated traits
#' * fraction_uncorr_traits: Fraction among the "nbr_traits" which are uncorrelated
#' * lambda: Integer, Pagel's lambda
#' * kappa: Integer, Pagel's kappa
#' * highCor: Correlation rate between the trait defined in "manTrait" and the simulated traits.
#' * manTrait: Index of the trait correlated with other traits.
#' @param param_tree list of integers needed for the simulation of the phylo object
#' @param dataframe data frame composed of 11 columns:
#' @param save path to save the data
#' @return list composed of a dataset with mixed data, a continuous data, a discrete data, the list of phylogenetic trees,
#' several list of the parameters used for the simulation and the dataframe used as input.
#' @export
data_simulator <- function(param_tree, dataframe, save = NULL){

  if(ncol(dataframe) != 11){
    stop("The number of columns should be equal to 11.")
  }

  #Rename columns of the data frame
  ################################
  newNames <-  c("nbr_traits", "class", "model", "states", "correlation",
                 "uncorr_traits", "fraction_uncorr_traits", "lambda", "kappa", "highCor", "manTrait")
  names(dataframe) <- newNames
  #in subclass there are, continuous, ordinal, interval(same quantity), nominal(no order).
  # uncorr_traits >= 0, 0<x<=1 fraction of uncovariant traits among the correlated group of traits, number of uncorrelated should be >=1.


  #Transform the columns (nbrs_traits as integer, uncorr_traits, fraction_uncorr_traits and lambda as numeric, and rest as character)
  #######################################################################
  dataframe[,c(2:3,5)] <- apply(dataframe[,c(2:3,5)], 2, as.character)
  dataframe$nbr_traits <- as.integer(dataframe$nbr_traits) #if the value is a float, converted in integer.
  dataframe[,c(4,6:11)] <- apply(dataframe[,c(4,6:11)], 2, as.numeric)


  #Check the data frame
  #####################

  #Homogenize the names
  dataframe$class[stringr::str_detect(dataframe$class, "^[cC]")] <- "continuous"
  dataframe$class[stringr::str_detect(dataframe$class, "^[oO]")] <- "ordinal"
  dataframe$class[stringr::str_detect(dataframe$class, "^[iI]")] <- "interval"
  dataframe$class[stringr::str_detect(dataframe$class, "^[nN]")] <- "non_eq_nominal"
  dataframe$model[stringr::str_detect(dataframe$model, "^[bB]")] <- "BM1"
  dataframe$model[stringr::str_detect(dataframe$model, "^[oO]")] <- "OU1"
  dataframe$model[stringr::str_detect(dataframe$model, "^[eE]")] <- "ER"
  dataframe$model[stringr::str_detect(dataframe$model, "^[sS]")] <- "SYM"
  dataframe$model[stringr::str_detect(dataframe$model, "^[aA]")] <- "ARD"
  dataframe$model[stringr::str_detect(dataframe$model, "^[nN]")] <- "Na"


  #Check if all the rows are filled correctly

  wrong <- which(dataframe$class == "continuous" &
                   (dataframe$model == "BM1" | dataframe$model == "OU1") & dataframe$states != 1)
  if(length(wrong) != 0){
    stop(paste("The line ", wrong, "is not filled correctly \n"))
  }

  wrong <- which(dataframe$class != "continuous" & dataframe$states <= 1)
  if(length(wrong) != 0){
    stop(paste("The line ", wrong, "is not filled correctly \n"))
  }

  wrong <- which(dataframe$lambda < 0 & dataframe$lambda > 1)
  if(length(wrong) != 0){
    stop(paste("The line ", wrong, "is not filled correctly \n"))
  }


  # Simulate phylogeny
  ####################
  birth = param_tree[[1]]
  death = param_tree[[2]]
  ntaxa = param_tree[[3]]

  # Simulating phylogenies fails sometimes. Try until we are successful
  Extant <- FALSE
  while (!Extant) {
    SimTree <- phytools::pbtree(b = birth, d = death, n = ntaxa, scale = 1, extant.only = TRUE)
    if(!is.null(SimTree)){
      Extant <- TRUE
    }
  }


  #Simulate by correlation
  ########################

  #Get the various group of correlated or not traits.
  correlation_values <- unique(dataframe$correlation)

  #Final matrix
  FinalData <- as.data.frame(matrix(0, nrow = length(SimTree$tip.label), ncol = 0))
  max_rate <- 0.5
  AlphasList <- list()
  ThetasList <- list()
  SigmasList <- list()
  TreeList <- list()
  NaParamList <- list()
  for(i in correlation_values){
    wrong <- which(sum(dataframe$correlation == i) == 1 & dataframe$correlation == i & dataframe$class != "continuous"
                   & dataframe$uncorr_traits != dataframe$nbr_traits &
                     dataframe$fraction_uncorr_traits != 0 & (dataframe$model != "OU1" | dataframe$model != "BM1"))
    if(length(wrong) != 0){
      stop(paste("The line ", wrong, "is not filled correctly \n"))
    }

    wrong <- which(dataframe$nbr_traits == 1 & dataframe$correlation == i & sum(dataframe$correlation == i) == 1 &
                     (dataframe$uncorr_traits != 1 | dataframe$fraction_uncorr_traits != 0))
    if(length(wrong) != 0){
      stop(paste("The line ", wrong, "is not filled correctly \n"))
    }

    #build a subset of traits being correlated together
    subdata <- subset(dataframe, dataframe$correlation == i)

    #extract rows in the dataframe corresponding to correlation group
    rows <- which(dataframe$correlation == i)
    #append row number in dataframe
    subdata <- cbind(subdata, rows)

    # #rescale phylogeny
    # lambdaCheck <- mean(subdata$lambda)
    # kappaCheck <- mean(subdata$kappa)
    #
    # if(lambdaCheck != 0 & kappaCheck != 0 & lambdaCheck != 1 & kappaCheck != 1){
    #   stop("lambda or kappa should be equal to 1")
    # }
    #
    # if(lambdaCheck != 1){
    #   subdataTree <- geiger::rescale(SimTree, "lambda", lambdaCheck)
    # }
    #
    # if(kappaCheck != 1){
    #   subdataTree <- geiger::rescale(SimTree, "kappa", kappaCheck)
    # }

    #check fraction of uncorrelated
    if(sum(subdata$fraction_uncorr_traits) != subdata$fraction_uncorr_traits[1] * nrow(subdata)){
      stop("The fraction of uncorrelated trait should be equal within the same group of correlated traits")
    }

    if(nrow(subdata) == 1 | (nrow(subdata) - sum(subdata$model == "Na") == 1)){

      if(nrow(subdata) == 1 && subdata$model == "Na"){
        stop("When model is Na, the traits should be correlated with another trait \n")
      }

      subdataOne <- subset(subdata, subdata$model != "Na")

      #rescale tree
      subdataTree <- rescaleTree(SimTree, subdataOne)

      Sigmas <- simSigma(subdataOne$nbr_traits, uncovTraits = subdataOne$uncorr_traits,
                         FracNocov = subdataOne$fraction_uncorr_traits)
      Thetas <- runif(subdataOne$nbr_traits, min = -10, max = 10)
      Alphas <- simAlpha(subdataOne$nbr_traits)
      #if(subdata$class == "continuous" | (subdata$class != "continuous" & subdata$uncorr_traits != subdata$nbr_traits)){
      if(subdataOne$class == "continuous" | (subdataOne$class != "continuous" &
                                             (subdataOne$model == "BM1" | subdataOne$model == "OU1"))){
        if(subdataOne$model == "BM1"){
          Alphas <- simAlpha(subdataOne$nbr_traits, alpha = 1e-5 * log(2))
        }


        #simulate independent continuous traits
        ContinuousData <- mvMORPH::mvSIM(tree = subdataTree,
                                         model = subdataOne$model,
                                         param = list(ntraits = subdataOne$nbr_traits,
                                                      theta = Thetas, #theta = ancestral states
                                                      sigma = Sigmas,
                                                      alpha = Alphas))


        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("F%s.%s/%s", 1:subdataOne$nbr_traits, i, subdataOne$row)

        if(any(subdataOne$class == "continuous")){
          FinalData <- cbind(FinalData, ContinuousData)
        }

        if(any(subdataOne$class != "continuous")){
          indexColumConvert <- 1:subdataOne$nbr_traits
          Nstates <- rep(subdataOne$states, subdataOne$nbr_traits)
          class <- rep(subdataOne$class, subdataOne$nbr_traits)
          ContinuousData <- scale(ContinuousData)

          DiscreteData <- ChangeContinuousTraitInDiscrete(ContinuousData, indexColumConvert, Nstates, class)

          #change columns names
          colnames(DiscreteData) <- sprintf("I%s.%s/%s", 1:subdataOne$nbr_traits, i, subdataOne$row)

          FinalData <- cbind(FinalData, DiscreteData)
        }
      }


      else{
        Ordinal <- FALSE
        if(subdataOne$class == "ordinal"){
          Ordinal <- TRUE
        }

        stateCheck <- 0
        #check if all the trait have the number of states requested
        while(stateCheck != subdataOne$nbr_traits){
          #simulate independent discrete traits
          DiscreteData <- simDiscreteTraits(subdataOne$nbr_traits,
                                            subdataOne$states,
                                            subdataOne$model,
                                            max_rate,
                                            subdataTree,
                                            equal = FALSE,
                                            Ordinal)

          stateCheck <- sum(apply(DiscreteData$tip_mat, 2,
                                  function(x){ifelse(length(unique(x)) < subdataOne$states, 0, 1)}))


          print(paste("ncol: ", stateCheck, "/", subdataOne$nbr_traits))

        }




        #change columns names
        DiscreteData$tip_mat <- as.data.frame(DiscreteData$tip_mat)
        colnames(DiscreteData$tip_mat) <- sprintf("I%s.%s/%s", seq(1:subdataOne$nbr_traits), i, subdataOne$row)


        FinalData <- cbind(FinalData, DiscreteData$tip_mat)

      }
    }

    if(nrow(subdata) > 1 & (nrow(subdata) - sum(subdata$model == "Na") > 1)){
      #check if the dataframe is correctly filled
      wrong <- which(subdata$model != "BM1" & subdata$model != "OU1")
      if(length(wrong) != 0){
        stop(paste("This line ", wrong, "is not filled correctly \n"))
      }

      ContiRows <- rep(subdata$rows[subdata$class == "continuous"],
                       subdata$nbr_traits[subdata$class == "continuous"]) #get the row index in the dataframe
      Ntraits <- sum(subdata$nbr_traits[ContiRows])
      model <- "BM1"
      indexBM <- NULL
      Thetas <- runif(Ntraits, min = -10, max = 10)
      Alphas <- simAlpha(Ntraits, alpha = 1e-5 * log(2)) #by default, alpha is for a Brownian motion model
      uncorrTraits <- sum(subdata$uncorr_traits[ContiRows])
      indexUncorrTraitsBM <- NULL
      uncorrTraitsBM <- sum(subdata$uncorr_traits[subdata$model == "BM1"])
      uncorrTraitsOU <- uncorrTraits - uncorrTraitsBM #number of uncorrelated traits simulated with OU
      #FracNocov <- sum(subdata$nbr_traits * subdata$fraction_uncorr_traits) / Ntraits
      Sigmas <- simSigma(Ntraits, uncovTraits = uncorrTraits, FracNocov = subdata$fraction_uncorr_traits[ContiRows])

      if("OU1" %in% subdata$model && "BM1" %in% subdata$model){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"

        # in case uncovTraits simulated with BM

        if(uncorrTraitsBM != 0){
          indexUncorrTraitsBM <- 1:uncorrTraitsBM
          indexBM <-c(indexUncorrTraitsBM, sample((1+uncorrTraits):Ntraits, (Ntraits - uncorrTraitsBM -
                                                                               sum(subdata$nbr_trait[subdata$model == "OU1"])), replace = FALSE))
        }
        if(uncorrTraitsBM == 0){
          #in case no uncovTraits in the BM model.
          indexBM <- sample((1 + uncorrTraitsOU):Ntraits,
                            Ntraits - sum(subdata$nbr_trait[subdata$model == "OU1"]), replace = FALSE)
        }
        # set 0 to diagonal for BM1
        diag(Alphas)[indexBM] <- 1e-5
      }

      if(sum(subdata$model == "OU1") == nrow(subdata)){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"
      }


      #rescale tree
      subdataTree <- rescaleTree(SimTree, subdata)

      #Simulate data
      ContinuousData <- mvMORPH::mvSIM(tree = subdataTree,
                                       model = model,
                                       param = list(ntraits = Ntraits,
                                                    theta = Thetas, #theta = ancestral states
                                                    alpha = Alphas))


      #If correlated discrete traits from continuous traits
      discreteSubdata <- subset(subdata, (subdata$class != "continuous" &
                                            (subdata$model == "BM1" | subdata$model == "OU1")))
      if(nrow(discreteSubdata) >= 1){
        Nstates <- rep(discreteSubdata$states, discreteSubdata$nbr_traits)
        class <- rep(discreteSubdata$class, discreteSubdata$nbr_traits)

        indexColumConvert <- c()
        for(r in 1:nrow(discreteSubdata)){
          if(discreteSubdata$model[r] == "BM1" & discreteSubdata$uncorr_traits[r] != 0){
            indexColumConvert <- c(indexColumConvert, tail(indexUncorrTraitsBM, n = discreteSubdata$uncorr_traits[r]))
            indexUncorrTraitsBM <- indexUncorrTraitsBM[-tail(indexUncorrTraitsBM, n = discreteSubdata$uncorr_traits[r])] #delete the index selected above
          }
          if(discreteSubdata$model[r] == "OU1" & discreteSubdata$uncorr_traits[r] != 0){
            indexUncorrTraitsOU <- (1 + uncorrTraitsBM):uncorrTraits
            indexColumConvert <- c(indexColumConvert, tail(indexUncorrTraitsOU, n = discreteSubdata$uncorr_traits[r]))
            indexUncorrTraitsOU <- indexUncorrTraitsOU[-tail(indexUncorrTraitsOU, n = discreteSubdata$uncorr_traits[r])] #delete the index selected above
          }
        }


        if(sum(discreteSubdata$model == "OU1") == nrow(discreteSubdata)){
          BMTraitsRemained <- NULL
          OUTraitsRemained <- sample(Ntraits[(uncorrTraitsOU+1):length(Ntraits)],
                                     sum(discreteSubdata$nbr_traits[discreteSubdata$model == "OU1"]) - uncorrTraitsOU)
        }

        if(sum(discreteSubdata$model == "BM1") == nrow(discreteSubdata)){
          BMTraitsRemained <- sample(Ntraits[(uncorrTraitsBM + 1):length(Ntraits)],
                                     sum(discreteSubdata$nbr_traits[discreteSubdata$model == "BM1"]) - uncorrTraitsBM)
          OUTraitsRemained <- NULL
        }

        indexOU <- NULL
        if("OU1" %in% discreteSubdata$model && "BM1" %in% discreteSubdata$model){
          BMTraitsRemained <- sample(indexBM[(uncorrTraitsBM + 1):length(indexBM)],
                                     sum(discreteSubdata$nbr_traits[discreteSubdata$model == "BM1"]) - uncorrTraitsBM)

          indexOU <- setdiff(1:Ntraits, indexBM)
          OUTraitsRemained <- sample(indexOU[(uncorrTraitsOU+1):length(indexOU)],
                                     sum(discreteSubdata$nbr_traits[discreteSubdata$model == "OU1"]) - uncorrTraitsOU)
        }


        if(!is.null(indexBM) & length(indexBM[(uncorrTraitsBM + 1):length(indexBM)]) == 1){
          BMTraitsRemained <- indexBM[(uncorrTraitsBM + 1):length(indexBM)]
        }


        if(!is.null(indexOU) & length(indexOU[(uncorrTraitsOU+1):length(indexOU)]) == 1){
          OUTraitsRemained <- indexOU[(uncorrTraitsOU+1):length(indexOU)]
        }

        indexColumConvert <- sort(c(indexColumConvert, BMTraitsRemained, OUTraitsRemained))
        DiscreteData <- ChangeContinuousTraitInDiscrete(ContinuousData, indexColumConvert, Nstates, class)

        #change columns names
        colnames(DiscreteData) <- sprintf("F%s.%s/%s", seq(1:sum(subdata$nbr_traits)), i, ContiRows)
        discreteRows <- rep(discreteSubdata$rows, discreteSubdata$nbr_traits)
        colnames(DiscreteData)[indexColumConvert] <- sprintf("I%s.%s/%s",
                                                             seq(1:sum(discreteSubdata$nbr_traits)), i, discreteRows)
        FinalData <- cbind(FinalData, DiscreteData)
      }

      else{

        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("F%s.%s/%s", seq(1:sum(subdata$nbr_traits)), i, ContiRows)
        FinalData <- cbind(FinalData, ContinuousData)
      }
    } #close condition if subdata >2

    if(any("Na" %in% subdata$model)){

      NaRow <- which(subdata$model == "Na")
      NaParam <- c()
      for(mR in 1:length(NaRow)){

        if(subdata$states[which(subdata$model == "Na" & subdata$class != "continuous")] !=
           subdata$states[which(subdata$model != "Na")]){
          stop("The Na trait must have the same number of states that the correlated discrete trait.")
        }

        manTraitValue <- subdata$manTrait[NaRow[mR]]

        if(manTraitValue == 0){
          manTraitValue <- sample(ncol(FinalData), 1)
        }


        if(subdata$class[NaRow[mR]] != "continuous"){

          DiscreteData <- corDiscTraitsOneTrait(subdata$nbr_traits[NaRow[mR]],
                                                FinalData[, manTraitValue],
                                                subdata$highCor[NaRow[mR]])

          NaParam <- c(NaParam, DiscreteData$param)
          DiscreteData <- DiscreteData$data

          #change columns names
          colnames(DiscreteData) <- sprintf("I%s.%s/%s",
                                            seq(1:subdata$nbr_traits[NaRow[mR]]),
                                            i,
                                            subdata$row[NaRow[mR]])

          FinalData <- cbind(FinalData, DiscreteData)
        }


        if(subdata$class[NaRow[mR]] == "continuous"){

          ContinuousData <- corContiTraitsOneTrait(subdata$nbr_traits[NaRow[mR]],
                                                   FinalData[, manTraitValue],
                                                   subdata$highCor[NaRow[mR]])

          NaParam <- c(NaParam, ContinuousData$param)
          ContinuousData <- ContinuousData$data

          #change columns names
          colnames(ContinuousData) <- sprintf("F%s.%s/%s",
                                              seq(1:subdata$nbr_traits[NaRow[mR]]),
                                              i,
                                              subdata$row[NaRow[mR]])

          FinalData <- cbind(FinalData, ContinuousData)
        }
      }
    }

    #Save parameters in an own list
    AlphasList[[i]] <- Alphas
    ThetasList[[i]] <- Thetas
    SigmasList[[i]] <- Sigmas
    TreeList[[i]] <- SimTree

    if(any("Na" %in% subdata$model)){
      NaParamList[[i]] <- NaParam
    }

    row.names(FinalData) <- SimTree$tip.label

    if(length(grep("F.", colnames(FinalData))) > 0){
      #Standardize continuous traits mean = 0  and sd = 1
      ContiIndex <- grep("F.", colnames(FinalData))
      FinalData[, ContiIndex] <- scale(FinalData[, ContiIndex, drop = FALSE])
    }

    if(length(grep("I.", colnames(FinalData))) > 0){
      #convert the discrete columns in factors
      DiscreteIndex<- grep("I.", colnames(FinalData))
      FinalData[ ,DiscreteIndex] <- lapply(FinalData[ ,DiscreteIndex, drop = FALSE], factor)

      #add levels in factor if number of levels = 1
      for (c in 1:ncol(FinalData[ ,DiscreteIndex, drop = FALSE])){
        if(length(levels(FinalData[ ,DiscreteIndex, drop = FALSE][,c])) == 1){
          if(levels(FinalData[ ,DiscreteIndex, drop = FALSE][,c]) == "0"){
            FinalData[ ,DiscreteIndex][,c] <- factor(FinalData[ ,DiscreteIndex, drop = FALSE][,c], levels = c("0", "1"))
          }
          else{
            colName <- names(FinalData[ ,DiscreteIndex])[c]
            row <- as.numeric(str_extract(colName, "(?<=\\/)\\d+"))
            Nstates <- dataframe[row, "states"]
            FinalData[ ,DiscreteIndex][,c] <- factor(FinalData[ ,DiscreteIndex, drop = FALSE][,c],
                                                     levels = as.character(0:(Nstates-1)))
          }
        }
      }
    }
  }

  #Define list of object
  ######################
  Data <- list(FinalData = FinalData,
               AlphaMatrices = AlphasList,
               Thetas = ThetasList,
               SigmaMatrices = SigmasList,
               TreeList = TreeList[1],
               PhyloParam = param_tree,
               Na_param = NaParamList,
               Dataframe = dataframe)

  #Save data
  ##########
  if(!is.null(save)){
    save(Data, file = paste0(save, ".RData"))

  }

  return(Data)
}
