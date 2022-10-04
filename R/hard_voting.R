#' @title Hard voting
#'
#' @description This function aggregates the output of several imputed dataframes selecting the most recurring value for
#' each row.
#'
#' @usage hard_voting(listOfdataframe)
#'
#' @param listOfdataframe list containing dataframes
#' @return A dataframe
#' @export
hard_voting <- function(listOfdataframe){
  HVdata <- datasets[[1]]

  for(d in 1:ncol(HVdata)){

    columnstoHV <- data.frame(matrix(NA, nrow = nrow(HVdata), ncol = length(datasets)))
    picked <- sapply(datasets,`[`, d)

    for(col in 1:length(picked)){
      columnstoHV[,col] <- picked[[col]]
    }

    #categorical variable
    if(is.factor(picked[[1]])){
      HVdata[,d] <- as.factor(apply(columnstoHV, 1, function(x) names(which.max(table(x)))))
    }
    #continuous variable
    else{
      HVdata[,d] <- as.numeric(apply(columnstoHV, 1, mean))
    }
  }

  return(HVdata)
}
