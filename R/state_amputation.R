#'@title Amputation of values according the state and the continuous value
#'
#'@description This function amputes value according their states and the continuous
#'values of the correlated trait.
#'
#'@usage state_amputation(data, state, amount, level, columnNa, columnConti)
#'
#'@param data dataframe
#'@param state character: name of the state in which ampute data
#'@param amount numeric: amount of missing value within the state
#'@param level character: if top, ampute highest values. If low, remove lowest values. If rdn, remove randomly values
#'@param columnNa integer: index of the column in which ampute data
#'@param columnConti integer: index of the column correlated to columnNa
#'
#'@return A dataframe with missing value in the columnNa
#'@export
state_amputation <- function(data,
                             state,
                             amount,
                             level,
                             columnNa,
                             columnConti){

  index_state <- which(data[, columnNa, drop = FALSE] == state)
  subData <- data[which(data[index_state, columnNa] == state), ]

  if(level == "rdn"){
    amputedData <- TDIP::mcar_miss_meca(amount, subData, 1)
  }

  else if(level == "top"){
    indexNa <- order(subData[,columnConti], decreasing = TRUE)[1:(nrow(subData[,columnConti, drop = FALSE])*amount)]
    amputedData <- subData
    amputedData[indexNa, columnNa] <- NA
  }

  else if(level == "low"){
    indexNa <- order(subData[,columnConti], decreasing = FALSE)[1:(nrow(subData[,columnConti, drop = FALSE])*amount)]
    amputedData <- subData
    amputedData[indexNa, columnNa] <- NA
  }

  data[index_state, ] <- amputedData
  return(data)

}

