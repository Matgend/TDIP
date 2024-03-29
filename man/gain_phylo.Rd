% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gain_phylo.R
\name{gain_phylo}
\alias{gain_phylo}
\title{GAIN}
\usage{
gain_phylo(missingData, variance_fraction, Data, batch_size = round(ncol(missingData)*0.2), hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL)
}
\arguments{
\item{missingData}{data.frame of 1 or more columns containing NAs}

\item{variance_fraction}{total amount of minimum variance to be represented by the eigenvectors which correspond to the
phylogenetic inertia}

\item{tree}{phylo object}

\item{batch_size}{integer}

\item{hint_rate}{numerical}

\item{alpha}{numerical, hyperparameter}

\item{epochs}{integer, iterations}

\item{hint}{dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)}
}
\value{
A list of list containing in the "tab" imputedData the imputed Data and in the "tab" parametersGain, the discriminant loss values, the generative loss values, the MSE loss values and the iterations correspond to these values (important to generate a plot).
}
\description{
This function imputes missing data for continuous and discrete traits applying the Generative Adversarial #'
Imputation Networks (GAIN)
}
