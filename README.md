# TDIP

Trait Data Imputation with Phylogeny (TDIP) is package allowing to the imputation of missing values in a trait dataset with the help or not of a phylogenetic tree. By using this package, the user can simulate trait datasets, include missing values according to several missing mechanisms, and / or apply various imputation approaches. 

```{r setup}
library(TDIP)
library(naniar)
```

# Installation

TDIP uses R and python. All the package needed can be installed via R. 

1. Install TDIP directly from Github using devtools
```{r, eval = FALSE}
install.packages("devtools")
library(devtools)
devtools::install_github("Matgend/TDIP")
```

2. To use GAIN, `python >= 3.9` needs to be installed as well as the `tensorflow`, `numpy`, `scikit-learn` and `tqdm` packages. The R package `reticulate` provide all the commands needed.
```{r, eval = FALSE}

#install 
install.packages("reticulate")
library(reticulate)
install_miniconda()

#install packages:

#Conda environement
reticulate::conda_install("EnvironmentName",c("tensorflow = 2.9.1", "numpy = 1.22.3", "scikit-learn =  1.1.1", "tqdm = 4.64.1"))

#Pip
reticulate::py_install(c("tensorflow = 2.9.1", "numpy = 1.22.3", "scikit-learn =  1.1.1", "tqdm = 4.64.1"), pip = TRUE)
```

# Data

The data can be simulated or empirical. To simulate a trait dataset, the function `data_simulator()` can be used. The function takes 2 arguments:

1. a dataframe composed of these columns:
* *nbr_traits*: Number of traits simulated with specific parameters
* *class*: Type of traits, (continuous, non_eq_nominal(categorical) and ordinal)
* *model*: Evolutionary model (BM1, OU1, ARD, SYM, ER, Manual)
* *states*: Number of states for discrete traits, (if continuous set it to 1)
* *correlation*: Index corresponding to the group of simulated traits which are correlated or not to other traits
* *uncorr_traits*: Among the "nbr_traits", it's the number of uncorrelated traits
* *fraction_uncorr_traits*: Fraction among the "nbr_traits" which are uncorrelated
* *lambda*: Pagel's lambda
* *kappa*: Pagel's kappa
* *highCor*: Correlation rate between the trait defined in "manTrait" and the simulated traits.
* *manTrait*: Index of the trait correlated with the trait(s) simulated according to the manual model.

2. A list containing the birth rate, the death rate and the number of taxa

The function `data_simulator()` returns a list containing the simulated trait data, the phylogenetic tree, all the parameters used for the simulation and the inputs.
```{r, results = "hide"}
#load the data present in the package
data(dataframe)
param_tree <- list(0.4, 0.1, 100)

simData <- data_simulator(param_tree = param_tree, 
                          dataframe = dataframe)

#print some rows of the simulated data
#head(simData$FinalData)
```


For an empirical dataset, you can use the function `data_preprocessing()` to convert all the character columns in factors and to create a list containing the dataset and the phylogenetic tree if provided.
```{r}

data(simulatedData)
data(tree)


simulatedData$Species <- rownames(simulatedData)

empData <- data_preprocessing(simulatedData, 
                              empTree = NULL, 
                              save = NULL)
#if no tree
empData <- data_preprocessing(simulatedData, 
                              empTree = tree, 
                              save = NULL)

#head(empData$FinalData)
```

# Generation of missing values

Missing values follow a missing mechanisms. In this package, the user can generate missing values according to 4 missing mechanisms:

1. Missing completely at random (MCAR)

2. Missing at random (MAR)

3. Missing not at random (MNAR)

4. PhyloNa
```{r}

missingRate <- 0.05

#MCAR
mcar_values <- mcar_miss_meca(missingRate = missingRate,
                              ds = simData$FinalData, cols_mis = 1:ncol(simData$FinalData))

#MAR
mar_values <- mar_miss_meca(missingRate = missingRate, 
                            ds = simData$FinalData, 
                            cols_mis = 1, 
                            cols_ctrl = 3)

#MNAR
mnar_values <- mnar_miss_meca(missingRate = missingRate, 
                              ds = simData$FinalData, 
                              cols_mis = colnames(simData$FinalData))

#PhyloNa
phyloNa_values <- phyloNa_miss_meca(missingRate = missingRate, 
                                    ds = simData$FinalData, 
                                    tree = simData$TreeList$`0`)
```

The function `na_insertion()` simulates missing values according to several missing mechanisms at once.
```{r, warning = FALSE}
missing_values <- na_insertion(missingRate = missingRate,
                               dataset = simData$FinalData, 
                               missingTraits = ncol(simData$FinalData), 
                               MARTraits = 1, 
                               MARctrlTraits = 3,
                               traitsNoNA = NULL, 
                               tree = simData$TreeList$`0`, 
                               save = NULL)


#visualize the amount of missing data
NaData <- missing_values$DataNaN

for(rdn in 1:length(NaData)){

  partition <- NaData[[rdn]][[1]]

  title <- names(partition)
  plot(vis_miss(partition), main = title)
}
```

# Missing values imputation

The imputation of missing values can be done by several imputation methods according several strategies. The methods available in the TDIP package are the following:

1. Phylogenetic imputation method
* Rphylopars `pi_continuous_traits_phylo()`
* corHMM `pi_categorical_traits_phylo()`

2. Machine learning
* missForest `missForest_phylo()`
* kNN `kNN_phylo()`
* predictive mean matching `mice_phylo()`

3. Deep learning
* GAIN `gain_phylo()`
* Polytomous regression `mice_phylo()`

Each of these methods can be applied according 3 strategies: 

1. no phylogenetic information

2. with phylogenetic information as eigenvectors

3. "2-step" which consists in using as input the fusion of a complete data obtained by using the phylogenetic imputation method to the incomplete data before the missing values imputation.

The function `missing_data_imputation()` provides a way to impute missing values through several imputation approaches and strategies. Each imputation method can be called separately too.
```{r, results = "hide"}
methods <- c("mice_phylo", "missForest_phylo", "kNN_phylo")
strategies <- "NP"
varfrac <- 0.95

#impute missing MCAR values
imputedData <- missing_data_imputation(ImputationApproachesNames = methods,
                                       data = missing_values$DataNaN$MCAR$`MCAR/13/0.05`, 
                                       tree = simData$TreeList$`0`, 
                                       strategies = strategies, 
                                       varfrac = varfrac, 
                                       save = NULL)
```

### Hard voting

Then, the function `hard_voting()` provides the possibility to apply an hard voting classifier using imputed datasets.
```{r}

#select imputed dataset in the imputedData object
namesMethods <- c("MICE", "MissForest", "KNN")
namesToselect <- paste(namesMethods, "NP", sep = "_")
namesImputedData <- names(imputedData)
datasets <- imputedData[which(namesImputedData %in% namesToselect)]

#hard voting
hv <- hard_voting(listOfdataframe = datasets)

```


# Error calculation

The error calculation is the root mean square error (RMSE) for continuous variables while for categorical variables is the proportion of falsely classified entries (PFC).
```{r}

#compute error for a data containing MCAR imputed by missForest
errors <- imputation_error(imputedData = imputedData$MissForest_NP,
                           trueData = simData$FinalData, 
                           missingData = missing_values$DataNaN$MCAR$`MCAR/13/0.05`, 
                           imputationApproachesName = methods[2], 
                           dataframe = simData$Dataframe)
```



