#' @title Eigenvectors calculations
#'
#' @description This function calculates the eigenvectors of a phylo object
#'
#' @usage get_eigenvec(tree, variance_fraction)
#'
#' @param tree phylogenetic tree of class "phylo"
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @return a data frame in which each column represents an eigenvector
#' @export
get_eigenvec <- function(tree, variance_fraction){

  decomp <- PVR::PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
  eigvec <- as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors

  egval <- decomp@Eigen$values #extract eigenvalues
  eigPerc <- egval/(sum(egval)) #calculate % of variance
  eigPercCum <- t(cumsum(eigPerc)) #cumulated variance

  #eigenvectors representing more than X% variance
  numEigen <- sum(eigPercCum < variance_fraction)

  eigOK <- eigvec[,1:numEigen, drop = T]

  if(numEigen == 0){
    print(paste("The variance_fraction should at leat be equal to ",
                eigPercCum[1]))
    eigOK <- eigvec[1]
  }
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  #Eigenvectors generated in object 'eigenTobind'
  #rename eigenTobind species column so it matches trait dataset species column
  eigOK <- as.data.frame(eigOK)
  names(eigOK)[1] <- "c1"
  row.names(eigOK) <- decomp@phylo$tip.label

  return(eigOK)
}
