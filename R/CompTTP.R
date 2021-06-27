CompTTP <- function(patdata, cwm=matrix(c(0, 0.1, 0.25, 0.5, 1, 10,     #patdata and not patdat
                                         0, 0.2, 0.5, 1, 2, 10,
                                         0, 0.2, 0.4, 1, NA, NA), byrow=TRUE, nrow=3)) {
  # Compute normalized Total Toxicity Profile (TTP).
  #
  # Args:
  #   patdat: patient toxicity data
  #   cwm: clinical weight matrix.
  #
  # Returns:
  #   A list of raw TTP, normalized TTP and an indicator of DLT.
  
  # convert toxicity data into matrix
  #toxdata <- data.matrix(patdata[ , c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5")]) - 1 
  toxdata <- matrix(as.numeric(as.matrix(patdata[ , c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5")])),nrow=dim(cwm)[1]) 
  # matrix of element prodcut
  epmat <- toxdata * cwm
  
  # compute TTP
  colmax <- apply(epmat, 1, max, na.rm=TRUE)
  ttp <- sqrt( t(colmax) %*% colmax )
  # whether there is DLT,
  if (max(colmax) > 1 ) dlt = 1
  else dlt = 0
  # compute maximum TTP
  uplim <- apply(cwm, 1, max, na.rm=TRUE)
  maxttp <- sqrt( t(uplim) %*% uplim )
  # normalized TTP
  nttp <- ttp / (maxttp + 0.0001)
  
  result <- list(ttp, nttp, dlt)
  names(result) <- c("TTP", "nTTP", "DLT")
  return(result) 
}
