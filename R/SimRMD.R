
SimRMD <- function(seed=2014, strDose=1, chSize=3, trlSize=36, numTrials=1000, sdose=1:6, MaxCycle=6, 
                        tox.target=0.28, control, iter=10000, burnin=4000, thin=1, chains=1, pathout="./", tox.matrix,
                        wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0  , 0   , 0.5, 1  ), 
                                        byrow = T, ncol = 5), toxmax = 2.5) {

# SimRMD <- function(seed=2441, strDose=1, chSize=3, trlSize=36, numTrials=1000, sdose=1:6, MaxCycle=6,
                        # tox.target=0.28, control, iter=10000, burnin=4000, thin=1, chains=1, pathout="./",
                        # sc, trend = 0.1) {						
						
# SimRMD <- function(strDose = 1, chSize = 3, trlSize = 36, numTrials =
                 # 1000, sdose = 1:6, MaxCycle = 6, sd.gamma = 5,
                 # sd.epsilon = 1, sc = probaT, formula = nTTP ~ dose +
                 # cycle, loss.fun = "default", tox.target = 0.28, trend
                 # = 0.1, seed = 2441, ...)

						
  # Simulate clinical trials performing repeated measures design (RMD).
  #
  # Args:
  #   seed: random seed for reproducing simulation results.
  #   strDose: starting dose level.
  #   chSize: cohort size of patient accrual.
  #   dlSize: dose level size limitation.
  #   trlSize: sample size of a clinical trial.
  #   numTrials: number of replicated clinical trials.
  #   iMax: maximum dose level.
  #   toxtype: toxicity types
  #   alpha: a vector of intercept with monotonic increasing order
  #   beta: slope for standardized dose
  #   sigma0: standard deviation of the random intercept
  #   sdose: a vector of standardized dose
  #   gamma: cycle effect
  #   cwm: clinical weight matrix.
  #   formula: formula in the linear mixed-effect model
  #   loss.fun: loss function for computing Bayes risk
  #   tox.target: the target toxicity/nTTP score.
  #   
  #   MaxCycle: define the maximum numbers of cycle.
  #
  # Returns:
  #   Simulation metrics for evaluating operating characteristics
  #     pctAlc: percentage of allocation.
  #     pctSlc: percentage of selection.
  #     actSize: actual sample size.
  #     numDLTs: number of DLTs.
  #     runTime: running time of the simulation.
options(warn=-1)
sd.gamma <- 5; sd.epsilon = 1;
formula <- nTTP ~ dose + cycle;
loss.fun <- "default";
set.seed(seed);  
  
if (loss.fun == 'default'){
	loss.fun <- function(data, nTTP, tox.target, dose, cycle){
	  # first test whether "cycle" is in the formula or not;
	  #idx.cycle <- match("cycle", colnames(fit$data$X.other))
	  idx.cycle <- match("cycle", all.vars(delete.response(terms(data$formula))))
	  if (is.na(idx.cycle)){
		idx.row <- (data$X.dose == dose)
	  }else{  
		idx.row <- (data$X.dose == dose) & (data$X.other[,idx.cycle] == cycle)
	  }
	  #print('nTTP');print(length(nTTP));print(nTTP);print(idx.row);
	  abs(nTTP[idx.row] - tox.target[1])
	}
}  
  
 #debug purpose
  #record patData
  rec_patData=NULL;
  rec_doseA=NULL;
  #masterData stores all the simulated y, dosage, cycle
  masterData=NULL;
  #debug  
  
  set.seed(seed)
  
  hyper <- list(control=control,iter=iter,burnin=burnin,thin=thin,chains=chains,pathout=pathout)
  #hyper <- list(control=control, iter=5000, burnin=500, thin=1, chains=1,pathout="res/")
  control <- hyper$control
  iter <- hyper$iter
  burnin <- hyper$burnin
  thin <- hyper$thin
  chains <- hyper$chains
  pathout <- hyper$pathout
  

  timeStart <- proc.time()  # set starting time
  
  
  # initialize storage matrices
  alcMat <- matrix(data=0, nrow=max(sdose)+1, ncol=numTrials)  # allocation
  recMat <- matrix(data=0, nrow=max(sdose)+1, ncol=numTrials)  # recommendation
  dropMat_DLT <- matrix(data=0, nrow=MaxCycle, ncol=numTrials)  # drop out by DLT event
  dropMat_RND <- matrix(data=0, nrow=MaxCycle, ncol=numTrials)  # drop out by RND event
  dropMat_all <- matrix(data=0, nrow=MaxCycle, ncol=numTrials)  # drop out by RND event
  nTTP_rec <- matrix(data=0, nrow=1, ncol=numTrials)  # nTTP at MTD
  nTTP_cycle_rec <- matrix(data=0, nrow=MaxCycle, ncol=numTrials)  # nTTP at MTD
  DLT_rec <- matrix(data=0, nrow=1, ncol=numTrials)  # DLT at MTD
  DLT_cycle_rec <- matrix(data=0, nrow=MaxCycle, ncol=numTrials)  # DLT at MTD
  obsDLT <- rep(0, length=numTrials)  # observed DLT
  actSize <- rep(0, numTrials)  # actual sample size
  doseTrace_AllTrials <- NULL;
  
  #reformat toxicity matrix
  dose.name <- NULL;
  for (d in sdose){dose.name <- c(dose.name, paste('Dose Level ',d,sep=''))}
  grade.name <- NULL;
  for (g in 1:5){grade.name <- c(grade.name, paste('Grade ',g-1,sep=''))}
  new.tox.matrix <- vector('list', MaxCycle);
  toxtype = c("Neurological", "Renal", "Hematological")
  for (c in 1:MaxCycle){
    new.tox.matrix[[c]] <- array(0, dim=c(length(toxtype), 5, length(sdose)),  dimnames = list(toxtype,grade.name,dose.name));
    for (d in sdose){
      this <- tox.matrix[d,c,,];
      new.tox.matrix[[c]][,,d] <- this;
    }
  }
  tox.matrix <- new.tox.matrix;
  
  # calculate the probability matrix for all dosage and cycles
  
  pb <- txtProgressBar(0,numTrials);
  for (k in 1:numTrials) {
    setTxtProgressBar(pb,k);
    currTime <- 0
    doseA <- strDose[1] # initial dosage
    cmPat <- chSize # initial sample size
    alcMat[doseA+1, k] <- chSize # initial dose allocation
    patData <- NULL
    doseTrace <- NULL
	masterData <- NULL;
    
    while (cmPat < trlSize)      
    {
      # simulate inter-arrival time for the coming cohort
      itarTime <- 1
      currTime <- currTime + itarTime
      
      ## assume 10 patients, each has 2 cycles
      n=3; K=MaxCycle;
      # design matrix W for the random effect
	  W <- NULL;
	  for (temp in 1:K){W <- rbind(W,diag(n))};
      dose.pat = rep(sdose[doseA], n)
      xx=expand.grid(dose.pat, 1:K)
      X=as.matrix(cbind(rep(1, n*K), xx))
      colnames(X) = c("intercept", "dose", "cycle")
	  #print('xx');print(dim(xx));print(xx);print(X);
	  #print('X');print(X);print('X_');
      # generate nTTP values
      #beta = beta;
      gamma = rnorm(n, sd=sd.gamma);
	  
	  #generate the nTTP based on the dose & cycle in X,
	  #for all the patients of one dose, all the cycle data are generated
	  #print('W');print(dim(W));
      #y = X %*% beta + W %*% gamma + rnorm(n*K,sd=sd.epsilon)
	  #generate the nTTP using Monia data
	  y = NULL;dlt=NULL;
	  for (i in 1:dim(X)[1]){
	    #temp=GenScn_nTTP(sc, dose=X[i,2], count.cycle=X[i,3],trend=trend);
		temp=GenTTP(tox.matrix, wm = wm, toxmax = toxmax, dose=X[i,2], count.cycle=X[i,3]);
		y=c(y,temp$nttp)
		dlt=c(dlt,temp$dlt);
		#print('test1');print(temp$nttp);print(temp$dlt);
	  }
	  #y = y + rnorm(n*K,sd=sd.epsilon);
	  #print('y');print(y);print(dim(W%*%gamma));
	  #masterData stores all the simulated y, 1, dosage, cycle, subject

	  masterData = rbind(masterData,cbind(y,X,rep(cmPat/chSize,n*K),1:n,dlt,runif(length(y))));
	  #print('masterData');print(masterData);


	#write the data from masterData to patData as input for the dosage estimation
	#writing depends on the current cohort, simulated cohort and simulated cycle
	#maximum writable cycle = current cohort - simulated cohort + 1
	#if simulated cycle <= maximum writable cycle, then write the cycle into dosage estimation input patData
	  current_cohort = cmPat/chSize
	  index = 1:dim(masterData)[1]
	  patData <- NULL
	#browser through all the data in the masterData
	  DLT_patientlist <- NULL;
	  drop_patientlist <- NULL;
      for (i in index){
		simulated_cohort = masterData[i,5];
		simulated_cycle = as.numeric(masterData[i,4]);
		simulated_subj = masterData[i,6];
		#ADD criteria:
		#simulated data are not available for patient with DLT
		if (paste(simulated_cohort,simulated_subj,sep='_') %in% DLT_patientlist){next;}
		if (paste(simulated_cohort,simulated_subj,sep='_') %in% drop_patientlist){next;}
		if (simulated_cycle <= ((current_cohort - simulated_cohort) + 1)){
			subID <- paste("cohort", simulated_cohort, "subj", simulated_subj, sep="")
			# indvData <- AcrPat(subID=subID, toxtype=c("Anemia", "Diarrhea", "Fatigue"), 
                          # alpha=c(1.5, 2.2, 3.0, 4.0, 5.0), beta=-1, sdose=c(1,2,3,4,5), sigma0=0.5,
                          # cdl=doseA, gamma=0.01, cycle=1) 
			# DLT <- any( apply(indvData[,6:10], 1, which.max) >= 4 )
            patData <- rbind(patData, 
                data.frame(uniqueID=subID, cohort=simulated_cohort, subj=simulated_subj, dose=masterData[i,3], cycle=simulated_cycle, nTTP=masterData[i,1], TTP=masterData[i,1], DLT=masterData[i,7]))
        # when generating data using LMM, nTTP = TTP = y
		}
	  	#record DLT patients
        if (masterData[i,7]==1){
			DLT_patientlist <- c(DLT_patientlist,paste(simulated_cohort,simulated_subj,sep='_'));
			if ((trlSize - chSize) == cmPat){
				dropMat_DLT[simulated_cycle,k] <- dropMat_DLT[simulated_cycle,k] + 1;
				dropMat_all[simulated_cycle,k] <- dropMat_all[simulated_cycle,k] + 1;
			}
		}
		#decide whether the patient will be dropped in the next cycle
		#random dropout criteria
		#10% dropout probability in the 2nd cycle, 20% in the 3rd cycle, 30% in the 4th cycle and more
		#dropout_prob <- simulated_cycle * 0.05;
		#if (simulated_cycle >= 3) {dropout_prob <- 0.15}
		#random dropout criteria 2		
		#disable random dropout probability in all cycles
		dropout_prob <- -1;
		if (masterData[i,8]<=dropout_prob){
			drop_patientlist <- c(drop_patientlist,paste(simulated_cohort,simulated_subj,sep='_'));
			if ((trlSize - chSize) == cmPat){
				dropMat_RND[simulated_cycle,k] <- dropMat_RND[simulated_cycle,k] + 1;
				if (!(paste(simulated_cohort,simulated_subj,sep='_') %in% DLT_patientlist)){
					dropMat_all[simulated_cycle,k] <- dropMat_all[simulated_cycle,k] + 1;
				}
			}
		}
	}
	#print('masterData');print(masterData);
	#print('patData');print(patData);
	#debug
	rec_patData <- rbind(rec_patData,patData,rep(-999,dim(patData)[2]));
	#debug

	
      doseTrace <- c(doseTrace, doseA)
      
      
      # making dose suggestion for the next cohort
	  #print('patData1');print(patData);
      nxtDose <- RunRMD2(formula=formula, data=patData, control=control, iter=iter, burnin=burnin, thin=thin, chains=chains, 
                        pathout=pathout, sdose=sdose, tox.target=tox.target, numTrials=numTrials, MaxCycle=MaxCycle, thisTrial=k)[2]

	  #update 10/14/2014
      #choke the dosage increase to 1
	  #print('nxtDose');print(nxtDose);print(nxtDose[1]);
	  if (as.numeric(nxtDose[1])>(doseA+1)){
		doseA<-doseA+1;
	  }else{
		doseA<-as.numeric(nxtDose[1]);
	  }

	  if (cmPat == chSize){
		if (sum(patData$DLT)==0){
			doseA <- 2;dose_flag=0;
		}else{
			doseA <- 1;dose_flag=1;			
		}
	  }
	  if ((cmPat >= chSize*2)&(dose_flag==1)){
		if (sum(patData$DLT[(length(patData$DLT)-chSize+1):length(patData$DLT)])==0){
			doseA <- 2;dose_flag=0;
		}else{
			doseA <- 1;dose_flag=1;			
		}
	  }	  
	  
      #doseA <- as.numeric(nxtDose[1])
	  rec_doseA <- c(rec_doseA, doseA);
		  
      cmPat <- cmPat + chSize # cmPat for the next cohort (loop)
      alcMat[doseA+1, k] <- alcMat[doseA+1, k] + chSize # aclMat for the next cohort (loop)
      
      # early stopping for over-toxicity 
      if (doseA == 0 ) {
        break;
        #warning("early stopping for over-toxicity!")
      }     
      
      # early stopping for exceeding upper limit for a dose level
      #if (max(alcMat[, k]) >= dlSize) {
      #  break;
      #}   
      
    }
    
    cmPat <- cmPat - chSize # cmPat for the next cohort (loop)
    alcMat[doseA+1, k] <- alcMat[doseA+1, k] - chSize # aclMat for the next cohort (loop)
    
    #plot(as.factor(1:(cmPat/chSize)), doseTrace, type='b', xlab="cohort", ylab="dose", main=paste("Trial", k))
	
	#save information for plotting all the dose Trace in one pdf file
    doseTrace_AllTrials[[k]]=doseTrace;
	
    #alcMat[doseA+1, doseB+1, k] <- alcMat[doseA+1, doseB+1, k] / cmPat
    recMat[doseA+1, k] <- recMat[doseA+1, k] + 1
    
	#write the nTTP and DLT of the 1st cycle of the MTD
	this_index <- (patData$dose == doseA) & (patData$cycle == 1);
	nTTP_rec[k] <- mean(patData$nTTP[this_index]);
	DLT_rec[k] <- mean(patData$DLT[this_index]);
	#record the nTTP and DLT of the each cycle
	for (this_cycle in 1:MaxCycle){
		this_index <- (patData$cycle == this_cycle);
		nTTP_cycle_rec[this_cycle,k] <- mean(patData$nTTP[this_index]);
		DLT_cycle_rec[this_cycle,k] <- sum(patData$DLT[this_index]);
	}
    obsDLT[k] <- sum(patData[,"DLT"])
    actSize[k] <- cmPat
    #try(write.table(recMat, file=paste(pathout, "/", k, "/recMat.txt", sep="")), silent=T)
	#try(write.table(alcMat, file=paste(pathout, "/", k, "/alcMat.txt", sep="")), silent=T)
    # write.table(dropMat_DLT, file=paste(pathout, "/", k, "/dropMat_DLT.txt", sep=""))   
	# write.table(dropMat_RND, file=paste(pathout, "/", k, "/dropMat_RND.txt", sep=""))
	# write.table(dropMat_all, file=paste(pathout, "/", k, "/dropMat_all.txt", sep=""))
	# write.table(nTTP_rec, file=paste(pathout, "/", k, "/nTTP_rec.txt", sep=""))
	# write.table(DLT_rec, file=paste(pathout, "/", k, "/DLT_rec.txt", sep=""))
	# write.table(nTTP_cycle_rec, file=paste(pathout, "/", k, "/nTTP_cycle_rec.txt", sep=""))
	# write.table(DLT_cycle_rec, file=paste(pathout, "/", k, "/DLT_cycle_rec.txt", sep=""))	
  }


  #debug
  #try(write.table(rec_patData,file=paste(pathout,'/patData.txt',sep='')), silent=T);
  #debug

  #debug
  #try(write.table(alcMat,file=paste(pathout,'/alcMat.txt',sep='')), silent=T);
  #debug

  #try(write.table(recMat,file=paste(pathout,'/recMat.txt',sep='')), silent=T);

  # write.table(dropMat_DLT,file=paste(pathout,'/dropMat_DLT.txt',sep=''));
  # write.table(dropMat_RND,file=paste(pathout,'/dropMat_RND.txt',sep=''));
  # write.table(dropMat_all,file=paste(pathout,'/dropMat_all.txt',sep=''));
  # write.table(nTTP_rec, file=paste(pathout, "/nTTP_rec.txt", sep=""))
  # write.table(DLT_rec, file=paste(pathout, "/DLT_rec.txt", sep=""))  
  # write.table(nTTP_cycle_rec, file=paste(pathout, "/nTTP_cycle_rec.txt", sep=""))
  # write.table(DLT_cycle_rec, file=paste(pathout, "/DLT_cycle_rec.txt", sep=""))	  
  
  #plot dose Trace
  #pdf(file=paste(pathout,'/doseTrace.pdf',sep=''));
  for (k in 1:numTrials) {
    plot(as.factor(1:length(doseTrace_AllTrials[[k]])), doseTrace_AllTrials[[k]], type='b', xlab="cohort", ylab="dose", main=paste("Trial", k),ylim=c(0,5))
  }
  dev.off();
  
  #timeLapse <- round((proc.time() - timeStart) / 60, 2)
  
  results <- list(alcMat/colSums(alcMat),
                  apply(recMat, 1, sum)/numTrials, 
                  median(obsDLT), 
                  median(actSize))
  names(results) <- c("Allocation", "Selection", 
                      "MedianDLT", "MedianSize"
                      )
  cat('\n');
  cat(paste('Operating characteristics based on ',numTrials,' simulations:\n',sep=''));cat('\n');
  cat(paste('Sample size ',trlSize,'\n',sep=''));cat('\n');
  output_matrix <- rbind(matrix(results$Allocation[2:(length(sdose)+1)],nrow=1), matrix(results$Selection[2:(length(sdose)+1)],nrow=1));
  dose.name <- NULL;
  for (d in sdose){dose.name <- c(dose.name, paste('Dose ',d,sep=''))}
  cycle.name <- NULL;
  for (c in 1:MaxCycle){cycle.name <- c(cycle.name, paste('mnTTP.',c,sep=''))}
  #print(output_matrix);print(dose.name)
  colnames(output_matrix) <- dose.name;
  rownames(output_matrix) <- c('Allocation %', 'Recommendation %');
  cat('$op.table\n');
  print(round(output_matrix,3));cat('\n');
  op.table <- output_matrix;
  #print out true nTTP
               nTTP.array <- function(wm) {
                        nTTP <- array(NA, c(rep(5, nrow(wm))))
                        capacity <- 5^(nrow(wm))
                        for(tt in 1 : capacity) {
                                index <- vec2array(tt, dim = c(rep(5, nrow(wm))))
                                nTTP[tt] <- sqrt(sum(wm[t(rbind(1 : nrow(wm), index))]^2)) / toxmax
                        }
                        return(nTTP)
                }
                nTTP.all <- nTTP.array(wm)

               pDLT <- function(proba, wm){    
                        DLT_array <- array(NA, c(rep(5, nrow(wm))))
                        capacity <- 5^(nrow(wm))
                        for(tt in 1 : capacity) {
                                index <- vec2array(tt, dim = c(rep(5, nrow(wm))))
                                DLT_array[tt] <- as.numeric(max(wm[t(rbind(1 : nrow(wm), index))]) >= 1)
                        }
                        return(sum(proba * DLT_array))
                }
                
                ####compute mean nTTP given probability array######
                mnTTP <- function(proba, nTTP.all){    
                        return(sum(proba * nTTP.all))
                }
                nTTP_prob <- function(nTTP.prob) {
                nTTProb <- array(NA, c(rep(5, nrow(nTTP.prob))))
                        capacity <- 5^(nrow(nTTP.prob))
                        for(tt in 1 : capacity) {
                                index <- vec2array(tt, dim = c(rep(5, nrow(nTTP.prob))))
                                nTTProb[tt] <- prod(nTTP.prob[t(rbind(1 : nrow(nTTP.prob), index))])
                        }
                        return(nTTProb)
                }
                #proba <- outer(outer(nTTP.prob[1, ], nTTP.prob[2, ]), nTTP.prob[3, ])
				sc.mat <- matrix(0,nrow=MaxCycle, ncol=max(sdose)); pdlt <- matrix(0,nrow=MaxCycle, ncol=max(sdose));
                for(i in 1:MaxCycle){   ##loop through cycle
                        for(j in 1:max(sdose)){
                                nTTP.prob <- tox.matrix[[i]][, ,j]
                                proba <- nTTP_prob(nTTP.prob)
								nTTP.all <- nTTP.array(wm)
                                sc.mat[i, j] <- round(mnTTP(proba, nTTP.all), 4)
                                pdlt[i, j] <- round(pDLT(proba, wm), 4)
                        }
                }
                sc <- rbind(sc.mat, pdlt[1, ])
                colnames(sc) <- dose.name
                rownames(sc) <- c(cycle.name,"pDLT")
				cat('$sc\n');
				sc <- round(sc, 3);
                print(sc);cat('\n');
				
				
  results <- list(alcMat/colSums(alcMat),
                  apply(recMat, 1, sum)/numTrials, 
                  median(obsDLT), median(actSize),
				  op.table, sc)
  names(results) <- c("Allocation", "Selection", 
                      "MedianDLT", "MedianSize", "op.table", "sc")
  return(results)
  
}
