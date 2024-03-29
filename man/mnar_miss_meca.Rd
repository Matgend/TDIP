% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnar_miss_meca.R
\name{mnar_miss_meca}
\alias{mnar_miss_meca}
\title{Impute NaNs in a simulated data according the MNAR mechanism}
\usage{
mnar_miss_meca(missingRate, ds, cols_mis)
}
\arguments{
\item{missingRate}{numerical vector corresponding to the rate of missing value to introduce in the data}

\item{ds}{dataframe in which NA should be created}

\item{cols_mis}{vector of index or columns name where missing value will be created}
}
\value{
A dataset with NA following a pattern of MNAR.
}
\description{
Simulate NaNs in a complete data composed of discrete and continuous traits which are correlated
or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MNAR. It used the function
delete_MNAR from the package missMethods.
}
