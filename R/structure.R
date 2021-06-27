#library(Formula)

set.data <- function(formula, data)
{
  ## remove the entry that has missing nTTP
  data <- data[!is.na(data$nTTP),]
  
  vals <- list(formula=formula)
  
  ## add the term "dose" if nTTP ~ cycle;
  if(is.na(match("dose", all.vars(formula))))
    formula <- update(formula, .~. + dose)
  
  vals$cohort.size <- nrow(data[(data$cohort == 1 & data$cycle == 1),])
  vals$n.cohort <- max(data$cohort)
  vals$n.cycle <- max(data$cycle)
    
  # number of subjects that has non-mising values
  vals$n.subj <- length(unique(data$uniqueID[!is.na(data$nTTP)]))
  
  
  # formula terms
  mf <- model.frame(formula, data)
  mt <- terms(mf)
  Y <- model.response(mf)
  X <- model.matrix(mt, mf)
  #print('data');print(data);
  vals$y <- Y
  
  # design matrix of random effects;
  pid.list <- unique(data$uniqueID) # unique patietn list, convert patient IDs to numeric orders;
  W <- matrix(0, length(vals$y), length(pid.list))
  for (i in 1:nrow(W)){
    data$pid[i] <- which(data[i, "uniqueID"]==pid.list)
    W[i, data$pid[i]] <- 1
    #print(data$pid[i])
  }
  vals$W <- W

  idx <- match("dose", all.vars(delete.response(mt))) # term index of p.var.bin in the right-hand side of the formula: 1;
  vals$idx.dose <- which(attr(X, "assign") == idx) # column index of p.var.bin in the design matrix X: 2
  vals$X.other <- as.matrix(X[, - (vals$idx.dose)]) # design matrix X without dose
  vals$P.other <- ncol(vals$X.other)
  vals$X.dose <- as.matrix(X[, vals$idx.dose])
  #print('X');print(idx);print('_');print(vals$idx.dose);print('_');print(X);
  #print('X.other');print(vals$X.other);print('P.other');print(vals$P.other);  
  vals
}

# parm <- function(prior = c("gamma", "invgamma", "normal"),
                 # monitor=TRUE, tuning=NULL, ...)
  #                 monitor=TRUE, center=FALSE, method=NULL, tuning=NULL, ...)

parm <- function(prior = c("gamma", "invgamma", "normal"), mean=0, var=100, shape=0.001, scale=0.001)
                 #monitor=TRUE, tuning=NULL, ...)
  #                 monitor=TRUE, center=FALSE, method=NULL, tuning=NULL, ...)
{
  monitor <- TRUE; tuning <- NULL;
  if (any(tuning <= 0)) stop("Tuning values must be positive")
  if (!is.logical(monitor)) stop("Monitor value must be a logical")
  #if (!is.logical(center)) stop("Center value must be a logical")
  
  retval <- list()
  retval$prior <- match.arg(prior)
  retval$monitor <- monitor
  #retval$center <- center
  #retval$method <- method
  retval$tuning <- tuning
  #hyper <- list(...)
  hyper <- list(mean=mean, var=var, shape=shape, scale=scale);
  
  switch(retval$prior,
         gamma = {
           if (!is.numeric(hyper$shape) || hyper$shape <= 0)
             stop("Gamma shape hyperparameter must be numeric > 0")
           if (!is.numeric(hyper$rate) || hyper$rate <= 0)
             stop("Gamma scale hyperparameter must be numeric > 0")
           retval$shape <- hyper$shape
           retval$rate <- hyper$rate
         },
         invgamma = {
           if (!is.numeric(hyper$shape) || hyper$shape <= 0)
             stop("Inverse gamma shape hyperparameter must be numeric > 0")
           if (!is.numeric(hyper$scale) || hyper$scale <= 0)
             stop("Inverse gamma scale hyperparameter must be numeric > 0")
           retval$shape <- hyper$shape
           retval$scale <- hyper$scale
         },
         normal = {
           if (!is.numeric(hyper$mean))
             stop("Normal mean hyperparameter must be numeric")
           if (!is.numeric(hyper$var))
             stop("Normal variance hyperparameter must be numeric")
           retval$mean <- hyper$mean
           retval$var <- hyper$var
           val <- try(chol(retval$var))
           # if (class(val) == "try-error")
           #   stop("Normal variance hyperparameter must be positive definite")
           retval$prec <- chol2inv(val)
           dim(retval$prec) <- dim(retval$var)
         }
  )
  
  structure(retval, class = "parm")
}

set.control <- function(control, sizes) {
  for(i in names(sizes)) {
    parm <- control[[i]]
    n <- sizes[[i]]
    switch(parm$prior,
           gamma = {
             parm$shape <- rep(parm$shape, length.out=n)
             parm$rate <- rep(parm$rate, length.out=n)
           },
           invgamma = {
             parm$shape <- rep(parm$shape, length.out=n)
             parm$scale <- rep(parm$scale, length.out=n)
           },
           normal = {
             parm$mean <- rep(parm$mean, length.out=n)
             if(is.matrix(parm$var) && any(dim(parm$var) != n))
               stop(paste0("Non-conformable covariance prior for '", i, "'"))
             else {
               parm$var <- diag(parm$var, n)
               parm$prec <- diag(parm$prec, n)
             }
           }
    )
    parm$n <- n
    control[[i]] <- parm
  }
  control
}

set.inits <- function(control)
{
  inits <- list()
  name=names(control)
  for(i in name) {
    parm <- control[[i]]
    switch(parm$prior,
           gamma = {
             vals <- rgamma(parm$n, 1, 1)
           },
           invgamma = {
             vals <- rgamma(parm$n, 1, 1)
           },
           normal = {
             vals <- rnorm(parm$n, 0, 1)
           },
		   #set the initial to true values for test
           normal = {
             vals <- rnorm(parm$n, 0, 0.001)
           },		   
           vals <- NULL
    )
    inits[[i]] <- vals
  }
  inits
}


is.parm <- function(x)
{
  class(x) == "parm"
}


Tox2nTTP <- function(tox = c(1,2,0), wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0  , 0   , 0.5, 1  ), 
                                        byrow = T, ncol = 5), toxmax = 2.5)
{
	res <- c(wm[,1][tox[1]],wm[,2][tox[2]],wm[,3][tox[3]]);
	return(sqrt(sum(res^2))/toxmax);
}
