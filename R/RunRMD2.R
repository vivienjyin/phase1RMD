 

RunRMD2 <- function(formula, data, control, iter, burnin, thin, chains, pathout, 
                   sdose, tox.target, numTrials, MaxCycle, thisTrial){
  
  fit <- BEst2(formula=formula, data=data, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains, pathout=pathout,thisTrial=thisTrial)
  if (thisTrial==1){
    fit2 <- BEst2(formula=formula, data=data, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains, pathout=pathout,thisTrial=thisTrial)
  }else{
    fit2 <- fit
  };
  #fit <- BEst2(formula=formula, data=patData, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains, pathout=pathout)
  try(summaryRMD(fit, fit2, pathout, numTrials, thisTrial), silent=T) ## TODO --> write the function to produce summary, effect size, traceplots, density plots;

  ## method 1: using existing observed patient covariates

  ###############################################################
  ## method 2: using simulated fake data for prediction
  
  # 1 patient a cohort
  # each patient has number of cycles that is equal to the maximum number of cycles in the observed data;
  
  data.rep <- NULL
  #for (i in 1:fit$data$n.cycle){
  for (i in 1:MaxCycle){
    for (j in 1:length(sdose)){
      subID <- paste("cohort", 1, "subj", j, sep="")
      data.rep <- rbind(data.rep, data.frame(subID, 1, j, sdose[j], i, 99999))
    }
  }
  #print('data.rep');print(data.rep);q(save='no');
  colnames(data.rep) <- c("uniqueID", "cohort", "subj", "dose", "cycle", "nTTP")
  fit$data <- set.data(formula, data.rep)
  #########################################################
  
  # compute the posteriors of nTTP
  nTTP.p <- NULL
  loss.DL <- NULL # loss at each dose level;
  #print('beta');print(fit$beta.other);
  #print('rmvnorm');
  #print(fit$gamma[1,]);print(fit$data$X.other);print(fit$data$X.dose);#print(fit$data$W);
  simu1=rnorm((iter-burnin)*length(sdose),mean=0, sd=rep(fit$s2.gamma[1:(iter-burnin),],each=length(sdose)))
  simu1=matrix(simu1,byrow=T,ncol=length(sdose));
  simu2=rnorm((iter-burnin)*dim(fit$data$X.dose)[1],mean=0, sd=rep(fit$s2.epsilon[1:(iter-burnin),],each=dim(fit$data$X.dose)[1]))
  simu2=matrix(simu2,byrow=T,ncol=dim(fit$data$X.dose)[1]);
  #print('fit$data$W');print(fit$data$W);print(simu1[1,]);
  for (i in 1:(iter-burnin)){
    mu.pred <- fit$data$X.dose %*% as.matrix(fit$beta.dose[i,]) + fit$data$X.other %*% as.matrix(fit$beta.other[i,])
	#print(fit$data$X.dose);print(dim(fit$data$X.other));print(fit$data$X.other);print(fit$beta.other[i,]);q(save='no');
    s2.pred <- fit$s2.gamma[i,] * tcrossprod(fit$data$W) + fit$s2.epsilon[i,] * diag(length(fit$data$y))
	#print('rmvnorm');
	#print(fit$gamma[i,]);print(fit$data$X.other);print(fit$data$X.dose);print(fit$data$W);
	#print(min(s2.pred));print(summary(c(s2.pred)));
	#print(summary(s2.pred[upper.tri(s2.pred)]));
	#print(summary(s2.pred[lower.tri(s2.pred)]));
	#save('mu.pred','s2.pred',file='debug.save.2.Rdata')
    #nTTP.pred <- rmvnorm(1, mean=mu.pred, sigma = s2.pred)
	
	nTTP.pred <- mu.pred + fit$data$W %*% matrix(simu1[i,],ncol=1);
	#print(rmvnorm(1,mean=rep(0,length(mu.pred)), sigma=diag(fit$s2.epsilon[i,],length(mu.pred))));
	nTTP.pred <- nTTP.pred + matrix(simu2[i,],ncol=1);
	#nTTP.pred <- mu.pred + fit$data$W %*% fit$gamma[i,];
	# define the loss function here;
  #loss.DC <- abs(nTTP.pred - (tox.target + fit$beta.other[i,2]*(fit$data$X.other[,2] - 1)))/ fit$data$X.other[,2] # loss at each dose & cycle
  loss.DC <- abs(nTTP.pred - (tox.target + 0*(fit$data$X.other[,2] - 1)))/ fit$data$X.other[,2] # loss at each dose & cycle
  #print('loss.DC');print(matrix(loss.DC, nrow=length(sdose), ncol=MaxCycle));
  #loss.D <- apply(matrix(loss.DC, nrow=length(sdose), ncol=MaxCycle), 1, sum)  # loss at each dose level;#average over cycle           
  loss.D <- matrix(loss.DC, nrow=length(sdose), ncol=MaxCycle)[,1]  # loss at the first dose level;#first cycle           
  #loss.D <- LossDose(fit, nTTP.pred, target.tox, beta.cycle=fit$beta.other[i,2], sDose, MaxCycle)
	
	#temp=matrix(0,nrow=dim(s2.pred)[1],ncol=dim(s2.pred)[2]);diag(temp)=diag(s2.pred);
	#nTTP.pred <- rmvnorm(1, mean=mu.pred, sigma = temp);
	#print('rmvnorm end');
    #nTTP.p <- cbind(nTTP.p, nTTP.pred)
	###!!!!!!new for cycle loss function###!!!!!!!
	#calculate cycle prediction
	#nTTP.pred <- rep(nTTP.pred,each=MaxCycle)+seq(0,MaxCycle-1)*fit$beta.other[i,2];
	###!!!!!!new for cycle loss function###!!!!!!!
  nTTP.p <- rbind(nTTP.p, nTTP.pred)
	loss.DL <- rbind(loss.DL, loss.D)
  }
  #print(cbind(nTTP.pred,fit$data$X.dose,fit$data$X.other,(tox.target + fit$beta.other[i,2]*(fit$data$X.other[,2] - 1))));
  #print('loss.DL');print(dim(loss.DL));
  #print('beta');print(mean(fit$beta.dose));
  #REVISE the target tox in multiple cycles with 
  #tox.target.new<-tox.target+seq(0,fit$data$n.cycle-1)*mean(fit$beta.other[1:(iter-burnin),2])
  #print('cycle');print(seq(0,MaxCycle-1)*fit$beta.other[i,2]);  
  #fit$nTTP.p <- apply(nTTP.p, 1, mean) # In method 1: nTTP.p is a 18 x 5000 matrix, fit$nTTP.p is a 18X1 vector;
  fit$nTTP.p <- apply(nTTP.p, 2, mean) # In method 1: nTTP.p is a 18 x 5000 matrix, fit$nTTP.p is a 18X1 vector;
  fit$Brisk <- apply(loss.DL, 2, mean)
  #print(fit$Brisk)
  #print('nTTP 0');print(matrix(fit$nTTP.p,nrow=length(tox.target)));
  #print(dim(nTTP.p))
  #print(length(fit$nTTP.p));q(save='no');
  # compute the expected loss (Bayes risk);
  # only look at cycle 1 first;
  #B.risk <- NULL
  #for (c in 1:1){
  #  for (d in sdose){
      #L <- loss.fun(data = fit$data, nTTP = fit$nTTP.p, tox.target = tox.target, dose = d, cycle = c)
	#  L <- loss.fun.cycle(data = fit$data, nTTP = fit$nTTP.p, tox.target = tox.target.new, dose = d, cycle = c, MaxCycle)
  #    B.risk <- rbind(B.risk, data.frame(mean(L), d, c))
  #  }
  #}
  #colnames(B.risk) <- c("BayesRisk", "dose", "cycle")
  #print('patData');print(data);save(data,file='patData.RData');
  #print('fit$data');print(summary(fit$beta.dose));print(summary(fit$beta.other[,1]));
  #print(summary(fit$beta.other[,2]));
  #print('B.risk');print(B.risk);
  
    # output the estimated dosage;
    min.br <- min(fit$Brisk)
    dose.est <- which.min(fit$Brisk)
    res <- c(min.br, dose.est)
    #print(res)
  #colnames(res) <- c("cycle", "min.Bayes.risk", "dose.suggestion")
  #print('dose rec');print(res);
  return(res)
}
    