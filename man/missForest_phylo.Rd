% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missForest_phylo.R
\name{missForest_phylo}
\alias{missForest_phylo}
\title{Non-parametric missing values imputation for mixed-type data by missForest(randomForest)}
\usage{
missForest_phylo(missingData, variance_fraction = 0, maxiter = 10, ntree = 100,
mtry = sqrt(ncol(missingData)), Data, hint = NULL)
}
\arguments{
\item{missingData}{data frame of 1 or more columns containing NAs}

\item{variance_fraction}{variance_fraction minimum variance (\%) explained by the eigenvectors}

\item{maxiter}{maximum number of iterations to be performed given the stopping criterion is not met beforehand.}

\item{ntree}{number of trees to grow in each forest.}

\item{mtry}{number of variables randomly sampled at each split. By default it's the square root of the number of
variables}

\item{tree}{phylo object}

\item{hint}{dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)}
}
\value{
A data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
}
\description{
This function imputes missing data for continuous and discrete traits applying the missForest approach
}
