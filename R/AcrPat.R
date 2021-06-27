AcrPat <- function(subID, toxtype, alpha, beta, sigma0, sdose, cdl, gamma, cycle) {
  # Simulate patient accrual with toxicity data specified by ToxProb()
  #
  # Args:
  #   subID: subject ID
  #   toxtype: toxicity types
  #   alpha: a vector of intercept with monotonic increasing order
  #   beta: slope for standardized dose
  #   sigma0: standard deviation of the random intercept
  #   sdose: a vector of standardized dose
  #   cdl: current dose level
  #   gamma: cycle effect
  #   cycle: the number of treatment cycle
  #
  # Returns:
  #   A data frame of patient toxicity data
  
  # number of toxicity types
  ntype <- length(toxtype)
  # initialize a matrix of patient toxicity data
  patTox <- matrix(NA, nrow=ntype, ncol=6)   ###ME: ncol can be variable, for example if we don't want to take into account the death because it requires another decision from clinicians

  celp <- ToxProb(toxtype, alpha, beta, sigma0, sdose, cdl, gamma, cycle)
  newcelp <- cbind(celp, 1 - apply(celp, 1, sum)) # adding grade 5  

  for (i in 1 : ntype) {
    patTox[i, ] <- stats::rmultinom(n=1, size=1, prob=newcelp[i, ]) #ME: is the same thing as sample(c(0:(5)),1, replace=TRUE,newcelp[i, ])?
  }  
  
  patData <- data.frame(cbind(rep(subID, ntype), rep(cdl, ntype ), rep(cycle, ntype), toxtype, patTox)) 
  colnames(patData) <- c("Pat ID", "Dose Level", "Cycle", "Type", "Grade 0", "Grade 1", "Grade 2", "Grade 3",
                         "Grade 4", "Grade 5")
  return(patData)
}