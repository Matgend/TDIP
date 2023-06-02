#'@title Amputation of values according the state and the continuous value
#'
#'@description Function to set NA for the same species in the random and top scenario
#'
#'@usage state_amputation_sameNA_top_rnd(data, state, amount, columnNa, columnConti)
#'
#'@param data dataframe
#'@param state character: name of the state in which ampute data
#'@param amount numeric: amount of missing value within the state
#'@param columnNa integer: index of the column in which ampute data
#'@param columnConti integer: index of the column correlated to columnNa
#'@return A dataframe with missing value in the columnNa
#'@export
state_amputation_sameNA_top_rnd <- function (data,
                                             state,
                                             amount,
                                             columnNa = 1,
                                             columnConti = 2) {
  DataRnd <- data
  DataTop <- data
  IdxC <- which(data[, 1] == state)
  NumNA <- round(amount * length(IdxC))
  SetNA <- sample(IdxC, NumNA)
  DataRnd[SetNA, columnNa] <- NA
  DataTop[SetNA, columnNa] <- NA
  Cont <- DataTop[IdxC, columnConti]
  Cont <- sort(Cont, decreasing = TRUE)
  DataTop[SetNA, columnConti] <- Cont[1:NumNA]
  DataTop[setdiff(IdxC, SetNA), columnConti] <- Cont[(NumNA + 1):length(Cont)]
  Out <- vector(mode = "list", length = 2)
  Out[[1]] <- DataRnd
  Out[[2]] <- DataTop
  return(Out)
}
