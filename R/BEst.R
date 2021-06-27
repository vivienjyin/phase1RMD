#04152015. change the M-H estimate of beta dose to gibbs sampling

######################################################
## model estimation using Bayesian methods
#####################################################

BEst <- function(formula, data, control, iter, burnin=0, thin=1, chains=1) {
  ptdata <- set.data(formula, data)
  X.dose <- ptdata$X.dose
  X.other <- ptdata$X.other
  W <- ptdata$W
  y <- ptdata$y
  
  control <- set.control(control, list(beta.dose=1, beta.other=ptdata$P.other, gamma=ptdata$n.subj, s2.gamma=1, s2.epsilon=1))
  
  ## Hyperparameters
  tau.beta <- control$beta.dose$tuning # increase the step of proposal, to jump to another value quicker;
  #u.beta <- control$beta.dose$shape; v.beta <- control$beta.dose$rate
  u.beta <- control$beta.dose$mean; v.beta <- matrix(1/control$beta.dose$var, 1,1)
  a.beta <- control$beta.other$mean; Bi.beta <- solve(control$beta.other$var)
  a.gamma <- control$s2.gamma$shape; b.gamma <- control$s2.gamma$scale
  a.epsilon <- control$s2.epsilon$shape; b.epsilon <- control$s2.epsilon$scale
  
  
  ## Model Parameters and Initial Values
  inits <- set.inits(control)
  beta.dose <- inits$beta.dose;
  if (beta.dose <0){beta.dose <-1};
  beta.other <- inits$beta.other;
  gamma <- inits$gamma
  s2.gamma <- inits$s2.gamma
  s2.epsilon <- inits$s2.epsilon
  
  res <- list()
  res$beta.dose <- matrix(NA, nrow = iter - burnin, ncol = 1)
  res$beta.other <- matrix(NA, nrow = iter - burnin, ncol = ptdata$P.other)
  res$s2.gamma <- matrix(NA, nrow = iter - burnin, ncol = 1)
  res$s2.epsilon <- matrix(NA, nrow = iter - burnin, ncol = 1)
  res$gamma <- matrix(NA, nrow = iter - burnin, ncol = ptdata$n.subj)
  
  ## if nTTP ~ dose + cycle;
  if(!is.na(match("dose", all.vars(ptdata$formula)))){
    
    accept <- 0
  #print('11');  
  beta.dose.rec <- NULL; mu.rec<-NULL;u.rec<-NULL; beta.cycle.rec<-NULL;
  for(i in 1:iter) {
      
    ## Beta for dose effect, using M-H algorithm
    # delta <- rnorm(1, 0, sd=tau.beta)
    # beta.dose.new <- exp(delta)*beta.dose
    #print('22');print(X.other);print('_');print(beta.other);
	#print('33');print(dim(X.dose));print(length(beta.dose.new));
	
	### This part is obsolete. do not use for M-H
    #ln.r <- dmvnorm(as.vector(y), mean=as.vector(X.dose %*% beta.dose.new + X.other %*% beta.other + W %*% gamma), sigma = s2.epsilon * diag(length(y)), log=TRUE) -
    #        dmvnorm(as.vector(y), mean=as.vector(X.dose %*% beta.dose + X.other %*% beta.other + W %*% gamma), sigma = s2.epsilon * diag(length(y)), log=TRUE) +
    #        dgamma(beta.dose.new, a.gamma, b.gamma, log=TRUE) - dgamma(beta.dose, a.gamma, b.gamma, log=TRUE) +
    #        delta
	### This part is obsolete: end
    
    # ln.r <- -1/(2*s2.epsilon) * (crossprod(y - W%*%gamma - X.other%*%beta.other - X.dose%*%beta.dose.new) - crossprod(y - W%*%gamma - X.other%*%beta.other - X.dose%*%beta.dose)) +
      # (u.beta-1) * (log(beta.dose.new) - log(beta.dose)) - v.beta * (beta.dose.new - beta.dose) +
      # delta
    
    # if (log(runif(1)) < ln.r) {
      # beta.dose <- beta.dose.new
      # accept <- accept + 1
    # }
    
    #print('12');
    ## The rest of the Beta, using Gibbs sampling
	
	##TEST purpose estimate beta 0 and beta cycle separately
	#print('s2');print(s2.epsilon);print(Bi.beta);print(a.beta);print(beta.other);
	# U <- chol(crossprod(X.other[,2]) / s2.epsilon + Bi.beta[2,2])
    # mu <- chol2inv(U) %*% (crossprod(X.other[,2], y - W %*% gamma - X.dose %*% beta.dose - X.other[,1] * beta.other[1]) / s2.epsilon + Bi.beta[2,2] * a.beta[2])
    # beta.other[2] <- mu + backsolve(U, rnorm(1))
	# #print('beta.other');print(mu);print(backsolve(U, rnorm(1)));print('end');

	# U <- chol(crossprod(X.other[,1]) / s2.epsilon + Bi.beta[1,1])
    # mu <- chol2inv(U) %*% (crossprod(X.other[,1], y - W %*% gamma - X.dose %*% beta.dose - X.other[,2] * beta.other[2]) / s2.epsilon + Bi.beta[1,1] * a.beta[1])
    # beta.other[1] <- mu + backsolve(U, rnorm(1))
	#print('beta.other[1]');print(beta.other);print(mean(gamma));print('end')

    #gibbs sampling of beta dose, constained to be > 0
	#print('X.dose');print(X.dose);print(s2.epsilon);print(v.beta);
    U <- chol(crossprod(X.dose) / s2.epsilon + v.beta)
    mu <- chol2inv(U) %*% (crossprod(X.dose, y - W %*% gamma - X.other %*% beta.other) / s2.epsilon + v.beta %*% u.beta)
    # print('beta.other');print(mu);print(backsolve(U, rnorm(ncol(X.other))));print('end');
	beta.dose.sample <- mu + rnorm(ncol(X.dose)*1000)/U
	temp=which(beta.dose.sample > 0)
	if (length(temp)>0){
		beta.dose=beta.dose.sample[temp[1]];
	}
    beta.dose.rec <- c(beta.dose.rec, beta.dose);
	beta.cycle.rec <- c(beta.cycle.rec, beta.other);
	mu.rec <- c(mu.rec, mu);
	u.rec <- c(u.rec, U);
	#print('debug X.dose');print(X.dose);print(beta.dose);  
	##Original code: estimate beta 0 and beta cycle together
    U <- chol(crossprod(X.other) / s2.epsilon + Bi.beta)
    mu <- chol2inv(U) %*% (crossprod(X.other, y - W %*% gamma - X.dose %*% beta.dose) / s2.epsilon + Bi.beta %*% a.beta)
	# print('beta.other');print(mu);print(backsolve(U, rnorm(ncol(X.other))));print('end');
    beta.other <- mu + backsolve(U, rnorm(ncol(X.other)))
	#print('debug X.dose end');
    #TEST estimate beta dose, beta 0 and beta cycle
	#print(v.beta);print(Bi.beta);print('end');
    #U <- chol(crossprod(cbind(X.dose,X.other)) / s2.epsilon + diag(c(v.beta,diag(Bi.beta))))
    #mu <- chol2inv(U) %*% (crossprod(cbind(X.dose,X.other), y - W %*% gamma ) / s2.epsilon + diag(c(v.beta,diag(Bi.beta))) %*% c(u.beta,a.beta))
	#print('beta.other');print(mu);print(backsolve(U, rnorm(ncol(X.other))));print('end');
    #temp <- mu + backsolve(U, rnorm(ncol(cbind(X.dose,X.other))))
	#beta.dose <- temp[1]; beta.other <- temp[2:3];
    
    ## Gamma
    U <- chol(crossprod(W) / s2.epsilon + diag(ncol(W)) / s2.gamma)
    mu <- chol2inv(U) %*% crossprod(W, y - X.dose %*% beta.dose - X.other %*% beta.other) / s2.epsilon
    gamma <- mu + backsolve(U, rnorm(ncol(W)))
    
    ## Random Effects Variance
    a <- 0.5 * length(gamma) + a.gamma
    b <- 0.5 * crossprod(gamma) + b.gamma
    s2.gamma <- 1 / rgamma(1, a, b)
	#print('s2.gamma');print(s2.gamma);
    
    ## Measurement Error Variance
    a <- 0.5 * length(y) + a.epsilon
    b <- 0.5 * crossprod(y - X.dose %*% beta.dose - X.other %*% beta.other - W %*% gamma) + b.epsilon
    s2.epsilon <- 1 / rgamma(1, a, b);
	#print('s2.epsilon');print(s2.epsilon);print(a);print(b);
	#print(y);print(X.dose);print(beta.dose);print(X.other);print(beta.other);print(b.epsilon);
    
    if (i > burnin) {
      res$beta.dose[i-burnin,] <- beta.dose
      res$beta.other[i-burnin,] <- beta.other
      res$s2.gamma[i-burnin,] <- s2.gamma
      res$s2.epsilon[i-burnin,] <- s2.epsilon
      res$gamma[i-burnin,] <- gamma
    }  
    
  }
  
  res$accept <- accept/iter
  }else{ 
  # if nTTP ~ cycle;
    
    for(i in 1:iter){
      ## The rest of the Beta, using Gibbs sampling
      U <- chol(crossprod(X.other) / s2.epsilon + Bi.beta)
      mu <- chol2inv(U) %*% (crossprod(X.other, y - W %*% gamma - X.dose %*% beta.dose) / s2.epsilon + Bi.beta %*% a.beta)
      beta.other <- mu + backsolve(U, rnorm(ncol(X.other)))
    
      ## Gamma
      U <- chol(crossprod(W) / s2.epsilon + diag(ncol(W)) / s2.gamma)
      mu <- chol2inv(U) %*% crossprod(W, y - X.dose %*% beta.dose - X.other %*% beta.other) / s2.epsilon
      gamma <- mu + backsolve(U, rnorm(ncol(W)))
    
    
      ## Random Effects Variance
      a <- 0.5 * length(gamma) + a.gamma
      b <- 0.5 * crossprod(gamma) + b.gamma
      s2.gamma <- 1 / rgamma(1, a, b)
      
      ## Measurement Error Variance
      a <- 0.5 * length(y) + a.epsilon
      b <- 0.5 * crossprod(y - X.dose %*% beta.dose - X.other %*% beta.other - W %*% gamma) + b.epsilon
      s2.epsilon <- 1 / rgamma(1, a, b)
    
      if (i > burnin) {
        res$beta.other[i-burnin,] <- beta.other
        res$s2.gamma[i-burnin,] <- s2.gamma
        res$s2.epsilon[i-burnin,] <- s2.epsilon
        res$gamma[i-burnin,] <- gamma
      }
    }
    
  }
  cmPat <- ptdata$n.subj
  # dir.create(paste(pathout, "/", thisTrial, sep=""));
  # write.table(beta.dose.rec, file=paste(pathout, "/", thisTrial, "/betadose",cmPat,".txt", sep=""))      
  # write.table(beta.cycle.rec, file=paste(pathout, "/", thisTrial, "/betacycle",cmPat,".txt", sep=""))      
  # write.table(mu.rec, file=paste(pathout, "/", thisTrial, "/mudose",cmPat,".txt", sep=""))      
  # write.table(u.rec, file=paste(pathout, "/", thisTrial, "/udose",cmPat,".txt", sep=""))      
  res$data <- ptdata
  return(res)
}

