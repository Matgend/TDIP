#' Simulated trait data set
#'
#' Data frame composed of 13 traits and 100 species. All the discrete traits are composed of 3 states. The trait simulated thought a brownian motion model or an MK model, are generated according an ultrametric tree with a strong phylogenetic signal.
#'
#' @format A data frame of 13 columns and 100 rows:
#' \describe{
#'   \item{I1.0/1}{Categorical trait simulated throught a MK model with a transition matrix obtained from an ARD rate model.}
#'   \item{I2.0/1, I3.0/1, I4.0/1}{Categorical traits simulated manually in a way that one state is highly correlated to a state in "I1.0/1".}
#'   \item{F1.0/3, F2.0/3, F3.0/3}{Continuous traits simulated throught a normal distribution and that are correlated at 80% with "I1.0/1".}
#'   \item{I1.0/4 I2.0/4 I3.0/4}{Categorical trait simulated throught a MK model with a transition matrix obtained from an ARD rate model. The 3 traits are independent.}
#'   \item{F1.1/2, F2.1/2, F3.1/2}{Continuous traits simulated throught a brownian motion model(BM). The 3 traits are independent.}
#' }
#' @source {Created in-house to serve as an example.}

#' @examples
#' data(simulatedData)
"simulatedData"


#' Simulated phylogenetic tree
#'
#' Ultrametric phylogenetic tree of class phylo simulated according a birth death process.
#'
#' @format A phylo object composed of 100 tips with a birth rate of 0.4 and a death rate of 0.1.
#' @source {Created in-house to serve as an example.}
#' @examples
#' data(tree)
"tree"


#' Input for the simulation of trait data
#'
#' Data frame which provides the structure of the data that the user wants to simulate.
#'
#' @format A data frame of 11 columns and 4 rows:
#' \describe{
#'   \item{nbr_traits}{Number of traits simulated with specific parameters.}
#'   \item{class}{Type of traits, (continuous, non_eq_nominal(categorical) and ordinal).}
#'   \item{model}{Evolutionary model (BM1, OU1, ARD, SYM, ER, Manual).}
#'   \item{states}{Number of states for discrete traits, (if continuous set it to 1).}
#'   \item{correlation}{Index corresponding to the group of simulated traits which are correlated or not to other traits.}
#'   \item{uncorr_traits}{Among the "nbr_traits", it's the number of uncorrelated traits.}
#'   \item{fraction_uncorr_traits}{Fraction among the "nbr_traits" which are uncorrelated.}
#'   \item{lambda}{Pagel's lambda.}
#'   \item{kappa}{Pagel's kappa.}
#'   \item{highCor}{Correlation rate between the trait defined in "manTrait" and the simulated traits.}
#'   \item{manTrait}{Index of the trait correlated with the trait(s) simulated according to the manual model.}
#' }

#' @source {Created in-house to serve as an example.}

#' @examples
#' data(dataframe)
"dataframe"
