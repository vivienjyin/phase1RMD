# prob <- GenScn_nTTP_probc('6', count.cycle = 2, trend = -0.1);
# prob_input <- NULL;
# prob_input$probaT1 <- prob[,,1];
# prob_input$probaT2 <- prob[,,2];
# prob_input$probaT3 <- prob[,,3];

gen_nTTP_DLT <- function(prob){
	prob_input <- NULL;
	prob_input$probaT1 <- prob[,,1];
	prob_input$probaT2 <- prob[,,2];
	prob_input$probaT3 <- prob[,,3];
	res <- nTTPbar(prob=prob_input);
	#the output is the nTTP and pDLT at each dose
	return(res);
}

GenCycle <- function(prob, count.cycle = 1, trend = 0.1){
  #print('prob');print(prob);
  cut = cutoff(prob)
  n.dose <- dim(prob$probaT1)[1]
  n.grade <- dim(prob$probaT1)[2]
  n.type = dim(cut)[3];
 
  #the nTTP of dose and cycle
  nTTP_res <- NULL;
  
  #res = nTTPbar(prob=prob)
  #res = nTTP_dist(prob=prob, dose=dose)
  
    a <- trend * (count.cycle - 1);
    proba.late = TP(cut, a) # shift mean to from 0 to a in the normal distribution to generate subsequent cycle tox prob;
    probc <- list(probaT1=proba.late[,,1], probaT2=proba.late[,,2], probaT3=proba.late[,,3])
	#resc = nTTP_dist(prob=probc, dose=dose)
 
  #print('probc');print(probc);
  return(proba.late)
}


GenAvgTTP <- function(prob){
	# prob_input <- NULL;
	# prob_input$probaT1 <- prob[,,1];
	# prob_input$probaT2 <- prob[,,2];
	# prob_input$probaT3 <- prob[,,3];
	#revision031917
	#res <- nTTPbar(prob=prob_input);
	res <- nTTPbar(prob=prob);
	#revision031917
	#the output is the nTTP and pDLT at each dose
	return(res);
}


GenTTP <- function(tox.matrix, wm, toxmax, dose=1, count.cycle = 1){
	prob <- tox.matrix[[count.cycle]][,,dose];
	#sample the toxicity grade for each toxicity type
	g <- NULL;
	for (t in 1:dim(prob)[1]){
		g <- c(g, which.max(rmultinom(1,1,prob[t,]))-1)
	}
	nttp <- Tox2nTTP(tox = g, wm = wm, toxmax = toxmax);
	dlt=0;dlt[(g[1] >= 3) | (g[2] >= 3) | (g[3] >= 4)]=1;
	res <- list(nttp = nttp, dlt = dlt);
	return(res);
}


ToxProb <- function(toxtype, alpha, beta, sigma0=0.5, sdose, cdl, gamma=0.01, cycle) {
  # Specify the toxicity probabilities for each type and grade at a given dose level using the proportional odds model
  #
  # Args:
  #   toxtype: a vector toxicity types in character string
  #   alpha: a vector of intercepts in monotonic increasing order
  #   beta: slope for standardized dose
  #   sigma0: standard deviation of the random intercept
  #   sdose: a vector of standardized doses
  #   cdl: current dose level
  #   gamma: slope for cycle
  #   cycle: the number of treatment cycle
  #
  # Returns:
  #   A matrix of toxicity probabilities (type times grade)
  
  if ( length(alpha) != 5 ) 
    stop("exactly five intercepts are needed for CTCAE grade 0--4!")
  else if ( min(diff(alpha, lag=1, differences=1)) < 0 )
    stop("intercepts must be in a monotonic increasing order!")
  
  # initialize storing matrix
  cump <- matrix(NA, nrow=length(toxtype), ncol=5)  # cumulative probability
  celp <- matrix(NA, nrow=length(toxtype), ncol=5)  # cell probability

  # a random effect from the same subject
  subj <- rnorm(1, mean=0 , sd=sigma0) 
  # an offset for various toxicity types
  offset <- runif(length(toxtype), min=-0.05, max=0.05)
  
  for (i in 1 : length(toxtype) ) {    
    # proportional odds model
    logitcp <- alpha + (beta + offset[i] ) * sdose[cdl] + subj + gamma * (cycle - 1)
    # cumulative probabilities
    cump[i, ] <- exp(logitcp) / (1 + exp(logitcp) ) 
    # cell probabilities
    celp[i, ] <- c(cump[i,1], diff(cump[i, ], lag=1, differences=1))
  }
  
  rownames(celp) <- toxtype
  colnames(celp) <- c('grade 0', 'grade 1', 'grade 2', 'grade 3', 'grade 4')
  return(celp)
}



# GenToxProb <- function(toxtype, intercept.alpha, coef.beta, sigma0=0.5, sdose = 1:length(intercept.alpha), cycle.gamma=0.01, cycle = 6) {
  # # Specify the toxicity probabilities for each type and grade at a given dose level using the proportional odds model
  # #
  # # Args:
  # #   toxtype: a vector toxicity types in character string
  # #   alpha: a vector of intercepts in monotonic increasing order
  # #   beta: slope for standardized dose
  # #   sigma0: standard deviation of the random intercept
  # #   sdose: a vector of standardized doses
  # #   cdl: current dose level
  # #   gamma: slope for cycle
  # #   cycle: the number of treatment cycle
  # #
  # # Returns:
  # #   A matrix of toxicity probabilities (type times grade)
  # dose.name <- NULL;
  # for (d in sdose){dose.name <- c(dose.name, paste('Dose Level ',d,sep=''))}
  # grade.name <- NULL;
  # for (g in 1:5){grade.name <- c(grade.name, paste('Grade ',g-1,sep=''))}
  # res <- vector('list', cycle);
  # for (c in 1:cycle){
    # res[[c]] <- array(0, dim=c(length(toxtype), 5, length(sdose)),  dimnames = list(toxtype,grade.name,dose.name));
    # for (d in sdose){
      # this <- ToxProb(toxtype = toxtype, alpha = intercept.alpha, beta = coef.beta, sigma0 = sigma0, sdose=sdose, cdl = d, gamma=cycle.gamma, cycle = c);
      # res[[c]][,,d] <- this;
    # }
  # }
  # return(res);
# }

#tox.matrix <- GenToxProb(toxtype = c("Renal", "Neuro", "Heme"), intercept.alpha = c(1.5, 2.2, 3, 4, 5), coef.beta = -1)
