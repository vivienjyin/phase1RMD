RunRMD <- function(data, control, trlSize = 36, tox.target=0.28, sdose = 1:6, strDose = 1, iter=10000, burnin=4000, thin=1, chains=1){

# RunRMD <- function(data, formula=nTTP ~ dose + cycle, tox.target=0.28, sdose=1:6, MaxCycle=6,
			# iter=10000, burnin=4000, thin=1, chains=1, seed = 2411, control=list(
			# beta.dose = parm("normal", mean = 0, var = 1000),
			# beta.other = parm("normal", mean = 0, var = 1000 ),
			# gamma = parm("normal", mean = 0, var = 100 ),
			# s2.gamma = parm("invgamma", shape = 0.001, scale = 0.001),
			# s2.epsilon = parm("invgamma", shape = 0.001, scale = 0.001))){


  formula <- nTTP ~ dose + cycle;
  strDose <- strDose; sdose <- sdose; 
  trlSize <- trlSize;
  MaxCycle <- 6;
  seed <- 2441;
  control <- list(
			beta.dose = parm("normal", mean = 0, var = 1000),
			beta.other = parm("normal", mean = 0, var = 1000 ),
			gamma = parm("normal", mean = 0, var = 100 ),
			s2.gamma = parm("invgamma", shape = 0.001, scale = 0.001),
			s2.epsilon = parm("invgamma", shape = 0.001, scale = 0.001));
			
  set.seed(seed);
  #data <- data$TTP;
  fit <- BEst(formula=formula, data=data, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains)
  # if (thisTrial==1){
    # fit2 <- BEst(formula=formula, data=data, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains)
  # }else{
    # fit2 <- fit
  # };
  fit2 <- fit;
  #fit <- BEst(formula=formula, data=patData, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains, pathout=pathout)
  #try(summaryRMD(fit, fit2, pathout, numTrials, thisTrial)) ## TODO --> write the function to produce summary, effect size, traceplots, density plots;
    
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
  nTTP.p <- cbind(nTTP.p, nTTP.pred)
	loss.DL <- cbind(loss.DL, loss.D)
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
  
  #summary the mean, median, SD, quantile
  q1fun <- function(x){quantile(x,0.025)};q2fun <- function(x){quantile(x,0.25)};
  q3fun <- function(x){quantile(x,0.75)};q4fun <- function(x){quantile(x,0.975)};
  mean.val <- apply(nTTP.p, 1, mean); median.val <- apply(nTTP.p, 1, median); SD.val <- apply(nTTP.p, 1, sd);
  q1 <- apply(nTTP.p, 1, q1fun); q2 <- apply(nTTP.p, 1, q2fun);
  q3 <- apply(nTTP.p, 1, q3fun); q4 <- apply(nTTP.p, 1, q4fun);
  #nTTP matrix MaxCycle * sdose
  #print('nTTP.pred');print(dim(nTTP.p));print(dim(loss.DL));
  #print('mean');print(dim(nTTP.p));print(mean.val);
  mean.matrix <- matrix(mean.val,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  median.matrix <- matrix(median.val,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  SD.matrix <- matrix(SD.val,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  q1.matrix <- matrix(q1,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  q2.matrix <- matrix(q2,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  q3.matrix <- matrix(q3,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  q4.matrix <- matrix(q4,byrow=T,nrow=MaxCycle,ncol=length(sdose));
  
  #print('fit$data');print(summary(fit$beta.dose));print(summary(fit$beta.other[,1]));
  #print(summary(fit$beta.other[,2]));
  #print('fit.nTTP');print(fit$nTTP.p);
  #print('nTTP.pred');print(nTTP.pred);print(nTTP[1,]);
    # output the estimated dosage;
	fit$Brisk <- abs(mean.matrix[1,] - tox.target);
    #print('fit.risk');print(fit$Brisk);
    min.br <- min(fit$Brisk)
    dose.est <- which.min(fit$Brisk)
    res <- c(min.br, dose.est)
  #colnames(res) <- c("cycle", "min.Bayes.risk", "dose.suggestion")
  #print('dose rec');print(res);
  dose.name <- NULL;
  for (d in sdose){dose.name <- c(dose.name, paste('Dose ',d,sep=''))}
  cycle.name <- NULL;
  for (d in 1:MaxCycle){cycle.name <- c(cycle.name, paste('toxpf',d,sep=''))}
  summary.matrix <- array(0, dim=c(MaxCycle,8,max(sdose)), dimnames = list(cycle.name,c('mean','sd','median','2.5%','25%','50%','75%','97.5%'),dose.name));
  for (i in 1:MaxCycle){
	  summary1.matrix <- rbind(mean.matrix[1,],SD.matrix[1,],median.matrix[1,]);
	  #print('summary1.matrix');print(summary1.matrix);
	  rownames(summary1.matrix) <- c('mean','sd','median');
	  colnames(summary1.matrix) <- paste0('Dose',sdose);
	  summary2.matrix <- rbind(q1.matrix[1,],q2.matrix[1,],median.matrix[1,],q3.matrix[1,],q4.matrix[1,]);
	  #print(summary2.matrix);
	  rownames(summary2.matrix) <- c('2.5%','25%','50%','75%','97.5%');
	  colnames(summary2.matrix) <- paste0('Dose',sdose);
	  summary.matrix[i,,] <- rbind(summary1.matrix, summary2.matrix);
  }
  #res <- list('nxtdose'=dose.est,'Estimate'=summary1.matrix, 'Quantiles'=summary2.matrix);
  cycle1.index <- seq(1, MaxCycle*length(sdose), MaxCycle);
  cycle1.nTTP <- nTTP.p[cycle1.index,];
  cat('\n');
  cat('Model: RMD with longitudinal toxicity\n');cat('\n');
  cat('Doses (skeleton):');cat('\n');
  print(sdose);
  cat('\n');
  cat(paste('The maximum sample size is: ',trlSize,sep=''));cat('\n');
  #cat(paste('The current enrolled number of patients are: ',max(patdata$cohort)*3,sep=''));cat('\n');
  #cat(paste('The current enrolled cohort is: ',max(patdata$cohort),sep=''));cat('\n');cat('\n');
  cat('Posterior estimates (mean) of toxicity:');cat('\n');
  colnames(mean.matrix) <- dose.name;rownames(mean.matrix) <- cycle.name;
  print(round(mean.matrix,3));cat('\n');
  cat(paste('Next recommended dose: ',dose.est,sep=''));cat('\n');
  res <- list('nxtdose'=dose.est,'tox.est'=round(summary.matrix,3), 'nTTP.p'=nTTP.p, 'sdose'=sdose, 'mean'=mean.matrix, 'MaxCycle'=MaxCycle);
  attr(res,'class') <- 'RunRMDVal'
  return(res)
}
