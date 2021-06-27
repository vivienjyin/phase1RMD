###################################################################################
###########Simulate Data from the Matrix and Do Stage 1-3##########################
###################################################################################
###model: yij = beta0 + beta1xi + beta2tj + ui + epiij#############################
######### ei = alpha0 + alpha1xi + alpha2xi^2 + gamma0ui + epi2i###################
######### epi ~ N(0, sigma_e); ui ~ N(0, sigma_u); epi2 ~ N(0, sigma_f)############
###################################################################################
#library(rjags)
#library(R2WinBUGS)
#library(BiocParallel)

# SimRMDEff <- function(seed=2441, tox_matrix, strDose = 1, chSize = 3, trlSize = 36, sdose = 1:6, MaxCycle = 6,  
                             # tox.target = 0.28, eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9)){

# # SimRMDEff <- function(seed=2441, tox_matrix, StrDose = 1, chSize = 3, trlSize = 36, sdose = 1:6, MaxCycle = 6,  
                             # # tox.target = 0.28, eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9),
                             # # eff.sd = 0.2,  p1 = 0.2, p2 = 0.2, 
                             # # ps1 = 0.2,  
                             # # proxy.thrd = 0.1, thrd1 = 0.28, thrd2 = 0.28, 
                             # # wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                            # # 0, 0.5, 0.75, 1  , 1.5, 
                                            # # 0, 0  , 0   , 0.5, 1  ), 
                                          # # byrow = T, ncol = 5), 
                             # # toxmax = 2.5, toxtype = NULL, intercept.alpha = NULL, 
                             # # coef.beta = NULL, cycle.gamma = NULL){							 
		
		# #revision031917
		# #dummy variable defination for winbugs parameters
		# beta_other <- 0;
		# beta_dose <- 0;
		# N1 <- 0;
		# inprod <- 0;
		# u <- 0;
		# N2 <- 0;
		# Nsub <- 0;
		# alpha <- 0;
		# gamma0 <- 0;
		# #revision031917

		
		# n.iter <- 3000;
		# eff.sd <- 0.2; p1 <- 0.2; p2 <- 0.2;
		# ps1 <- 0.2; proxy.thrd <- 0.1; thrd1 <- 0.28; thrd2 <- 0.28;
		# wm <- matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                            # 0, 0.5, 0.75, 1  , 1.5, 
                                            # 0, 0  , 0   , 0.5, 1  ), 
                                          # byrow = T, ncol = 5);
		# toxmax <- 2.5; toxtype = NULL; intercept.alpha = NULL;
		# coef.beta <- NULL; cycle.gamma <- NULL;
		
		# set.seed(seed)        
		# trialSize <- trlSize;
		# doses <- sdose;
		# cycles <- 1:MaxCycle;
		
        # listdo <- paste0("D", doses)
        # if(length(eff.structure) != length(doses))
                # stop("Check if you have specified the efficacy mean structure for the correct number of doses.")
        # if(nrow(wm) != 3)
                # stop("Make sure that you specified the clinical weight matrix for three types of toxicities.")
        
        # MaxCycle <- length(cycles)
        # Numdose <- length(doses)
        # flag.tox.matrix.null <- 0
        
        # if(is.null(tox_matrix)){
                # flag.tox.matrix.null <- 1
                # tox_matrix <- array(NA, dim = c(Numdose, MaxCycle, 3, 5))
                # if(length(toxtype) != 3)
                        # stop("Right now we are only considering three toxicity types, but we will relax this constraint in the future work.")
                # if(length(intercept.alpha) != 4){
                        # stop("Exactly four intercepts alpha are needed for grade 0--4 in simulation!")
                # } 
                # if(min(diff(intercept.alpha, lag = 1, differences = 1)) < 0){
                        # stop("Intercepts alpha for simulation must be in a monotonic increasing order!")
                # }
                # if(length(coef.beta) != 3){
                        # stop("Exactly three betas need to be specified for three types of toxicities!")
                # }
        # }
        
                
        
        # #############################################
        # #########Determin Scenarios##################
        # #############################################
        # w1 <- wm[1, ]
        # w2 <- wm[2, ]
        # w3 <- wm[3, ]
        
        # ####Array of the normalized scores########
        # nTTP.array <- function(w1,w2,w3){
                # nTTP <- array(NA,c(5,5,5))
                # for (t1 in 1:5){
                        # for (t2 in 1:5){
                                # for (t3 in 1:5){
                                        # nTTP[t1,t2,t3] <-  
                                                # sqrt(w1[t1]^2 + w2[t2]^2 + w3[t3]^2)/toxmax
                                # }
                        # }
                # }
                # return(nTTP)
        # }
        # nTTP.all <- nTTP.array(w1,w2,w3)
        
        # ####compute pdlt given probability array####
        # pDLT <- function(proba){    
                # return(sum(proba[c(4,5),      , ]) 
                       # + sum(proba[      ,c(4,5), ]) 
                       # + sum(proba[      ,      ,5])
                       # - sum(proba[c(4,5),,])*sum(proba[,c(4,5),])
                       # - sum(proba[c(4,5),,])*sum(proba[,,5])
                       # - sum(proba[,c(4,5),])*sum(proba[,,5])
                       # + sum(proba[c(4,5),,])*sum(proba[,c(4,5),])*sum(proba[,,5]))
        # }
        
        # ####compute mean nTTP given probability array######
        # mnTTP    <-function(proba){    
                # return(sum(proba * nTTP.all))
        # }
        
        # sc.mat <- matrix(NA, MaxCycle, Numdose)
        # pdlt <- matrix(NA, MaxCycle, Numdose)
        
        # for(i in 1:MaxCycle){   ##loop through cycle
                # for(j in 1:Numdose){ ##loop through dose
                        # if(flag.tox.matrix.null == 0){
                                # nTTP.prob <- tox_matrix[j, i, ,]
                        # }else{
                                # # initialize storing matrix
                                # cump <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cumulative probability
                                # celp <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cell probability
                                
                                # for(k in 1:length(toxtype)){    
                                        # # proportional odds model
                                        # logitcp <- intercept.alpha + coef.beta[k] * j +  cycle.gamma * (i - 1)
                                        # # cumulative probabilities
                                        # cump[k, 1:4] <- exp(logitcp) / (1 + exp(logitcp)) 
                                        # cump[k, 5] <- 1
                                        # # cell probabilities
                                        # celp[k, ] <- c(cump[k,1], diff(cump[k, ], lag = 1, differences = 1))
                                # }
                                
                                # nTTP.prob <- celp
                                # tox_matrix[j, i, ,] <- celp
                        # }
                        
                        # proba <- outer(outer(nTTP.prob[1, ], nTTP.prob[2, ]), nTTP.prob[3, ])
                        # sc.mat[i, j] <- round(mnTTP(proba), 4)
                        # pdlt[i, j] <- round(pDLT(proba), 4)
                # }
        # }
        
        # mEFF <-  eff.structure
        # sc <- rbind(sc.mat, mEFF, pdlt[1, ])
        # colnames(sc) <- listdo
        # rownames(sc) <- c("mnTTP.1st", "mnTTP.2nd", "mnTTP.3rd", 
                          # "mnTTP.4th", "mnTTP.5th", "mnTTP.6th",
                          # "mEFF", "pDLT")
        
        # ####################################################################################


        # k <- 1  ####Let us use k to denote the column of altMat####
        # doseA <- strDose ####Let us use doseA to denote the current dose assigned
        # cmPat <- chSize ####current patient size
        # masterData <- list()
        # rec_doseA <- NULL
        # stage1.allow <- NULL
        
        # altMat <- matrix(0, Numdose, k)
        # altMat[doseA, k] <- chSize
        
        # ###############################################################
        # ########################Stage 1################################
        # ###############################################################
        # #print('Stage 1')
        # while(cmPat <= trialSize/2){
                # rec_doseA <- c(rec_doseA, doseA)
                # PatData <- NULL
                # n <- chSize
                # K <- MaxCycle
                
                # ##############################################################
                # #############design matrix for toxicity#######################
                # ##############################################################
                # dose.pat <- rep(doseA, n)
                # xx <- expand.grid(dose.pat, 1:K)
                # X <- as.matrix(cbind(rep(1, n * K), xx))
                # X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
                # colnames(X) <- c("intercept", "dose", "cycle", "icycle")
                
                # ##############################################################
                # #######design matrix for efficacy#############################
                # ##############################################################
                # Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
                # colnames(Xe) <- c("intercept", "dose", "squared dose")
                
                # #####simulate nTTP and Efficacy for the cohort########
                # outcome <- apply(X, 1, function(a){
                        # tox.by.grade <- tox_matrix[a[2], a[3], ,]
                        # nttp.indices <- apply(tox.by.grade, 1, function(o){
                                # sample(1:5, 1, prob = o)
                        # })
                        # if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                                # y.dlt <- 1
                        # }else{
                                # y.dlt <- 0
                        # }
                        # y.nTTP <- nTTP.all[nttp.indices[1], 
                                           # nttp.indices[2],
                                           # nttp.indices[3]]
                        # return(c(y.dlt = y.dlt,
                                 # y = y.nTTP))
                # })
                # y <- outcome[2, ]
                # y.dlt <- outcome[1, ]
                
                # beta.mean <- eff.structure[Xe[1, 2]]
                # beta.sd <- eff.sd
                # beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
                # beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
                # e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
                
                
                # #####construct master Dataset and extract from this dataset for each interim#########
                # #####PatData stores the dataset used to estimate the model at each interim###########
                # temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n * K), 1:n))
                # temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
                # masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
                
                # temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
                # masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
                
                # current_cohort <- cmPat/chSize
                # PatData.index <- data.frame(cohort = 1:current_cohort, cycles = current_cohort:1)
                # PatData.index <- within(PatData.index, cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles))
                # for(i in 1:nrow(PatData.index)){
                        # PatData <- rbind(PatData, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                              # masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                      # c(1, 2, 3, 4, 5, 6, 9)])
                # }
                
                # #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt#######
                # dlt.drop <- sapply(by(PatData, PatData$subID, function(a){a[, 1]}), function(o){
                        # if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                # return(rep(0, length(o)))
                        # }else{
                                # o.temp <- rep(0, length(o))
                                # o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                # return(o.temp)
                        # }
                # })
                
                # dropout.dlt <- NULL
                
                # for(i in unique(PatData$subID)){
                        # dropout.dlt[PatData$subID == i] <- dlt.drop[[i]]
                # }
                
                # PatData <- PatData[dropout.dlt == 0, ]
                
                # ###################################################################################
                # ################Model Estimation given PatData#####################################
                # ################Stage 1 only considers the toxicity model##########################
                # ###################################################################################
                # nTTP <- PatData[, 2]
                # X_y <- as.matrix(PatData[, c("intercept", "dose", "icycle")])
                # n.subj <- length(unique(PatData[, "subID"]))
                # W_y <- matrix(0, nrow(PatData), n.subj)
                # W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$subID, function(a){
                        # which(a == unique(PatData$subID))
                # })))] <- 1
                
                # model.file <- function()
                # {       
                        # beta <- c(beta_other[1], beta_dose, beta_other[2])
                        # for(i in 1:N1){
                                # y[i] ~ dnorm(mu[i], tau_e)
                                # mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                        # }
                        # for(j in 1:N2){
                                # u[j] ~ dnorm(0, tau_u)
                        # }
                        # beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                        # beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                        # tau_e ~ dunif(p1_tau_e, p2_tau_e)
                        # tau_u ~ dunif(p1_tau_u, p2_tau_u)
                # }
                
                # mydata <- list(N1 = length(nTTP), N2 = n.subj, y = nTTP, X_y = X_y, 
                               # W_y = W_y, p1_beta_other = c(0, 0), 
                               # p2_beta_other = diag(rep(0.001, 2)),
                               # p1_beta_dose = 0, p2_beta_dose = 1000,
                               # p1_tau_e = 0, p2_tau_e = 1000, 
                               # p1_tau_u = 0, p2_tau_u = 1000)
                
                # path.model <- file.path(tempdir(), "model.file.txt")
                # R2WinBUGS::write.model(model.file, path.model)
                
                # inits.list <- list(list(beta_other = c(0.1, 0.1),
                                        # beta_dose = 0.1,
                                        # tau_e = 0.1,
                                        # tau_u = 0.1,
                                        # .RNG.seed = sample(1:1e+06, size = 1), 
                                        # .RNG.name = "base::Wichmann-Hill"))
                
                # jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                             # quiet = TRUE, inits = inits.list)
                
                # update(jagsobj, n.iter = n.iter, progress.bar = "none")
                # post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other"), 
                                                    # n.iter = n.iter, 
                                                    # progress.bar = "none")
                
                # if(cmPat == trialSize/2){
                        # ######################################################################
                        # #############define allowable doses###################################
                        # ######################################################################
                        # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                     # post.samples$beta_other[2,,1]))
                        # ####condition 1####
                        # prob1.doses <- sapply(doses, function(a){
                                # mean(apply(sim.betas, 2, function(o){
                                        # as.numeric(o[1] + o[2] * a <= thrd1)
                                # }))
                        # })
                        
                        # ####condition 2####
                        # prob2.doses <- sapply(doses, function(a){
                                # mean(apply(sim.betas, 2, function(o){
                                        # as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                                # }))
                        # })
                        
                        # allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                        # if(length(allow.doses) == 0){
                                # stage1.allow <- 0
                                # return(results = list(sc = sc, opt.dose = 0, 
                                                      # altMat = altMat, 
                                                      # masterData = masterData, 
                                                      # stage1.allow = stage1.allow, 
                                                      # dose_trace = rec_doseA, 
                                                      # cmPat = cmPat))
                        # }
                        # stage1.allow <- allow.doses
                        
                # }else{
                        # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1]))
                        # loss.doses <- sapply(doses, function(a){
                                # mean(apply(sim.betas, 2, function(o){
                                        # abs(o[1] + o[2] * a - tox.target)
                                # }))
                        # })
                        # nxtdose <- doses[which.min(loss.doses)]
                        
                        # if(as.numeric(nxtdose) > (doseA + 1)){
                                # doseA <- doseA + 1;
                        # }else{
                                # doseA <- as.numeric(nxtdose);
                        # }
                        
                        # if (cmPat == chSize){
                                # dlt <- with(masterData$toxicity, y.dlt[cycle == 1 & V7 == cmPat/chSize])
                                # if (sum(dlt) == 0){
                                        # doseA <- 2
                                        # dose_flag <- 0
                                # }else{
                                        # doseA <- 1
                                        # dose_flag <- 1			
                                # }
                        # }
                        # if ((cmPat >= chSize*2) & (dose_flag == 1)){
                                # dlt <- with(masterData$toxicity, y.dlt[cycle == 1 & V7 == cmPat/chSize])
                                # if (sum(dlt) == 0){
                                        # doseA <- 2
                                        # dose_flag <- 0
                                # }else{
                                        # doseA <- 1
                                        # dose_flag <- 1			
                                # }
                        # }
                        # if(cmPat >= chSize*3){
                                # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                             # post.samples$beta_other[2,,1]))
                                # ####condition 1####
                                # prob1.doses <- sapply(doses, function(a){
                                        # mean(apply(sim.betas, 2, function(o){
                                                # as.numeric(o[1] + o[2] * a <= thrd1)
                                        # }))
                                # })
                                
                                # ####condition 2####
                                # prob2.doses <- sapply(doses, function(a){
                                        # mean(apply(sim.betas, 2, function(o){
                                                # as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                                        # }))
                                # })
                                
                                # allow.doses <- which(prob1.doses >= ps1 & prob2.doses >= ps1)
                                # if(length(allow.doses) == 0){
                                        # stage1.allow <- 0
                                        # return(results = list(sc = sc, opt.dose = 0, 
                                                              # altMat = altMat, 
                                                              # masterData = masterData, 
                                                              # stage1.allow = stage1.allow, 
                                                              # dose_trace = rec_doseA, 
                                                              # cmPat = cmPat))
                                # }
                        # }
                        # altMat[doseA, k] <- altMat[doseA, k] + chSize
                # }
                # cmPat <- cmPat + chSize # cmPat for the next cohort (loop)
        # }
        # #print('Stage 2')
        # ###############################################################
        # ########################Stage 2################################
        # ###############################################################
        # #######let us wait until the efficacy data for stage 1 is available######
        # #########################################################################
        # cmPat <- cmPat - chSize
        
        
        # PatData <- list() 
        # PatData.index <- data.frame(cohort = current_cohort:1, cycles = 3:(3 + current_cohort - 1))
        # PatData.index <- within(PatData.index, cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles))
        # for(i in 1:nrow(PatData.index)){
                # PatData$toxicity <- rbind(PatData$toxicity, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                                        # masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                                # c(1, 2, 3, 4, 5, 6, 9)])
        # }
        
        # #####change the ordering of the subjects--being consistent######
        # temp.toxicity <- NULL
        # for(i in 1:current_cohort){
                # temp.toxicity <- rbind(temp.toxicity, PatData$toxicity[grepl(paste0("cohort",i,"subject"), PatData$toxicity$subID), ])
        # }
        
        # PatData$toxicity <- temp.toxicity
        # rm(temp.toxicity)
        
        # #####toxicity dropout########
        # #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
        # dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                # if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                        # return(rep(0, length(o)))
                # }else{
                        # o.temp <- rep(0, length(o))
                        # o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                        # return(o.temp)
                # }
        # })
        
        # dropout.dlt <- NULL
        
        # for(i in unique(PatData$toxicity$subID)){
                # dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
        # }
        
        # PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
        
        # #####efficacy dropout########
        # dropout.eff <- sapply(dlt.drop, function(a){
                # if(sum(a) == 0){
                        # return(0)
                # }else if(min(which(a == 1)) >= 4){
                        # return(0)
                # }else{
                        # return(1)
                # }
        # })
        
        # dropout.eff <- as.numeric(dropout.eff)
        # PatData$efficacy <- masterData$efficacy[dropout.eff == 0, ]
        
        
        # ############Model Estimation to randomize among allowable doses##########
        # nTTP <- PatData$toxicity[, 2]
        # X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
        # n.subj <- length(unique(PatData$toxicity[, "subID"]))
        # W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
        # W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                # which(a == unique(PatData$toxicity$subID))
        # })))] <- 1
        
        # EFF <- PatData$efficacy[, 1]
        # X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
        # keepeff.ind <- which(dropout.eff == 0)
        
        # model.file <- function()
        # {       
                # beta <- c(beta_other[1], beta_dose, beta_other[2])
                # for(i in 1:N1){
                        # y[i] ~ dnorm(mu[i], tau_e)
                        # mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                # }
                # for(k in 1:Nsub){
                        # u[k] ~ dnorm(0, tau_u)
                # }
                # for(j in 1:N2){
                        # e[j] ~ dnorm(mu_e[j], tau_f)
                        # mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                # }
                # beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                # beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                # alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                # gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                # tau_e ~ dunif(p1_tau_e, p2_tau_e)
                # tau_u ~ dunif(p1_tau_u, p2_tau_u)
                # tau_f ~ dunif(p1_tau_f, p2_tau_f)
        # }
        
        # mydata <- list(N1 = length(nTTP), N2 = length(EFF), y = nTTP, 
                       # Nsub = n.subj, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                       # W_y = W_y, X_e = X_e, p1_beta_other = c(0, 0), 
                       # p2_beta_other = diag(rep(0.001, 2)),
                       # p1_beta_dose = 0, p2_beta_dose = 1000,
                       # p1_alpha = c(0, 0, 0), p2_alpha = diag(rep(0.001, 3)), 
                       # p1_gamma0 = 0, p2_gamma0 = 0.001,
                       # p1_tau_e = 0, p2_tau_e = 1000, 
                       # p1_tau_u = 0, p2_tau_u = 1000, 
                       # p1_tau_f = 0, p2_tau_f = 1000)
        
        # path.model <- file.path(tempdir(), "model.file.txt")
        # R2WinBUGS::write.model(model.file, path.model)
        
        # inits.list <- list(list(beta_other = c(0.1, 0.1),
                                # beta_dose = 0.1,
                                # alpha = c(0.1, 0.1, 0.1),
                                # gamma0 = 0.1,
                                # tau_e = 0.1,
                                # tau_u = 0.1,
                                # tau_f = 0.1,
                                # .RNG.seed = sample(1:1e+06, size = 1), 
                                # .RNG.name = "base::Wichmann-Hill"))
        
        # jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                     # quiet = TRUE, inits = inits.list)
        
        # update(jagsobj, n.iter = n.iter, progress.bar = "none")
        # post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                       # "gamma0"), 
                                            # n.iter = n.iter, 
                                            # progress.bar = "none")
        
        # ######################################################################
        # #############update allowable doses###################################
        # ######################################################################
        # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                     # post.samples$beta_other[2,,1]))
        # ####condition 1####
        # prob1.doses <- sapply(doses, function(a){
                # mean(apply(sim.betas, 2, function(o){
                        # as.numeric(o[1] + o[2] * a <= thrd1)
                # }))
        # })
        
        # ####condition 2####
        # prob2.doses <- sapply(doses, function(a){
                # mean(apply(sim.betas, 2, function(o){
                        # as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                # }))
        # })
        
        # allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
        # if(length(allow.doses) == 0){
                # return(results = list(sc = sc, opt.dose = 0, 
                                      # altMat = altMat, 
                                      # masterData = masterData, 
                                      # stage1.allow = stage1.allow, 
                                      # dose_trace = rec_doseA, 
                                      # cmPat = cmPat))
        # }
        
        # sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
        # RAND.EFF <- sapply(allow.doses, function(a){
                # mean(apply(sim.alphas, 2, function(o){
                        # as.numeric(o[1] + o[2] * a + o[3] * a^2)
                # }))
        # })
        # RAND.EFF <- exp(RAND.EFF)/sum(exp(RAND.EFF))
        # nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)
        
        # ####a condition for untried higher dose level that is randomized#### rec_doseA[length(rec_doseA)]
        # if(nxtdose > max(rec_doseA) + 1 & !nxtdose %in% rec_doseA){
                # nxtdose <- max(rec_doseA) + 1
        # }
        
        
        # #################################################################
        # ########generate data for the new enrolled cohort in Stage 2#####
        # #################################################################
        # while(cmPat < trialSize - chSize){
                # doseA <- nxtdose
                # rec_doseA <- c(rec_doseA, doseA)
                # altMat[doseA, k] <- altMat[doseA, k] + chSize
                # cmPat <- cmPat + chSize
                # PatData <- list()
                # n <- chSize
                # K <- MaxCycle
                
                # #######################################################
                # #################design matrix for toxicity############
                # #######################################################
                # dose.pat <- rep(doseA, n)
                # xx <- expand.grid(dose.pat, 1:K)
                # X <- as.matrix(cbind(rep(1, n*K), xx))
                # X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
                # colnames(X) <- c("intercept", "dose", "cycle", "icycle")
                
                # #########################################################
                # #######design matrix for efficacy########################
                # #########################################################
                # Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
                # colnames(Xe) <- c("intercept", "dose", "squared dose")
                
                # #####simulate nTTP and Efficacy for the cohort########
                # outcome <- apply(X, 1, function(a){
                        # tox.by.grade <- tox_matrix[a[2], a[3], ,]
                        # nttp.indices <- apply(tox.by.grade, 1, function(o){
                                # sample(1:5, 1, prob = o)
                        # })
                        # if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                                # y.dlt <- 1
                        # }else{
                                # y.dlt <- 0
                        # }
                        # y.nTTP <- nTTP.all[nttp.indices[1], 
                                           # nttp.indices[2],
                                           # nttp.indices[3]]
                        # return(c(y.dlt = y.dlt,
                                 # y = y.nTTP))
                # })
                # y <- outcome[2, ]
                # y.dlt <- outcome[1, ]
                
                # beta.mean <- eff.structure[Xe[1, 2]]
                # beta.sd <- eff.sd
                # beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
                # beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
                # e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
                # #################################################################
                # #################################################################
                
                # #####add to the master Dataset and extract from this dataset for each interim#########
                # #####PatData stores the dataset used to estimate the model at each interim############
                # temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n*K), 1:n))
                # temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
                # masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
                
                # temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
                # masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
                
                # current_cohort <- cmPat/chSize
                # PatData.index <- data.frame(cohort = 1:current_cohort, cycles = current_cohort:1)
                # cycles.adj <- c(rep(2, trialSize/(2 * chSize)), rep(0, current_cohort - trialSize/(2 * chSize)))
                # PatData.index <- within(PatData.index, {cycles <- cycles + cycles.adj
                # cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles)})
                # for(i in 1:nrow(PatData.index)){
                        # PatData$toxicity <- rbind(PatData$toxicity, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                                                # masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                                        # c(1, 2, 3, 4, 5, 6, 9)])
                # }
                
                # #####toxicity dropout########
                # #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
                # dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                        # if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                # return(rep(0, length(o)))
                        # }else{
                                # o.temp <- rep(0, length(o))
                                # o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                # return(o.temp)
                        # }
                # })
                
                # dropout.dlt <- NULL
                # for(i in unique(PatData$toxicity$subID)){
                        # dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
                # }
                
                # PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
                # cohort.eff.index <- with(PatData.index, which(cycles >= 3))
                
                # #####efficacy dropout########
                # dropout.eff <- sapply(dlt.drop, function(a){
                        # if(sum(a) == 0){
                                # return(0)
                        # }else if(min(which(a == 1)) >= 4){
                                # return(0)
                        # }else{
                                # return(1)
                        # }
                # })
                
                # dropout.eff <- as.numeric(dropout.eff)
                
                # PatData$efficacy <- subset(masterData$efficacy, masterData$efficacy[, 5] %in% cohort.eff.index & dropout.eff == 0)
                
                
                
                # #############################################################################################################################
                # ######################Model Estimation to update toxicity information and randomization probs################################
                # #############################################################################################################################
                # nTTP <- PatData$toxicity[, 2]
                # X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
                # n.subj <- length(unique(PatData$toxicity[, "subID"]))
                # W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                # W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                        # which(a == unique(PatData$toxicity$subID))
                # })))] <- 1
                
                # EFF <- PatData$efficacy[, 1]
                # X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
                # keepeff.ind <- which(dropout.eff == 0 & masterData$efficacy[, 5] %in% cohort.eff.index)
                
                # model.file <- function()
                # {       
                        # beta <- c(beta_other[1], beta_dose, beta_other[2])
                        # for(k in 1:Nsub){
                                # u[k] ~ dnorm(0, tau_u)
                        # }
                        # for(i in 1:N1){
                                # y[i] ~ dnorm(mu[i], tau_e)
                                # mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                        # }
                        # for(j in 1:N2){
                                # e[j] ~ dnorm(mu_e[j], tau_f)
                                # mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                        # }
                        
                        # beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                        # beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                        # alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                        # gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                        # tau_e ~ dunif(p1_tau_e, p2_tau_e)
                        # tau_u ~ dunif(p1_tau_u, p2_tau_u)
                        # tau_f ~ dunif(p1_tau_f, p2_tau_f)
                # }
                
                # mydata <- list(N1 = length(nTTP), N2 = length(EFF), Nsub = n.subj, 
                               # y = nTTP, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                               # W_y = W_y, X_e = X_e, p1_beta_other = c(0, 0), 
                               # p2_beta_other = diag(rep(0.001, 2)),
                               # p1_beta_dose = 0, p2_beta_dose = 1000,
                               # p1_alpha = c(0, 0, 0), p2_alpha = diag(rep(0.001, 3)), 
                               # p1_gamma0 = 0, p2_gamma0 = 0.001,
                               # p1_tau_e = 0, p2_tau_e = 1000, 
                               # p1_tau_u = 0, p2_tau_u = 1000, 
                               # p1_tau_f = 0, p2_tau_f = 1000)
                
                # path.model <- file.path(tempdir(), "model.file.txt")
                # R2WinBUGS::write.model(model.file, path.model)
                
                # inits.list <- list(list(beta_other = c(0.1, 0.1),
                                        # beta_dose = 0.1,
                                        # alpha = c(0.1, 0.1, 0.1),
                                        # gamma0 = 0.1,
                                        # tau_e = 0.1,
                                        # tau_u = 0.1,
                                        # tau_f = 0.1,
                                        # .RNG.seed = sample(1:1e+06, size = 1), 
                                        # .RNG.name = "base::Wichmann-Hill"))
                
                # jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                             # quiet = TRUE, inits = inits.list)
                
                # update(jagsobj, n.iter = n.iter, progress.bar = "none")
                # post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                               # "gamma0"), 
                                                    # n.iter = n.iter, 
                                                    # progress.bar = "none")
                
                # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                             # post.samples$beta_other[2,,1]))
                
                # sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                
                # ############redefine allowable doses#############
                # ####condition 1####
                # prob1.doses <- sapply(doses, function(a){
                        # mean(apply(sim.betas, 2, function(o){
                                # as.numeric(o[1] + o[2] * a <= thrd1)
                        # }))
                # })
                
                # ####condition 2####
                # prob2.doses <- sapply(doses, function(a){
                        # mean(apply(sim.betas, 2, function(o){
                                # as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                        # }))
                # })
                
                # allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                # if(length(allow.doses) == 0){
                        # return(results = list(sc = sc, opt.dose = 0, 
                                              # altMat = altMat, 
                                              # masterData = masterData, 
                                              # stage1.allow = stage1.allow, 
                                              # dose_trace = rec_doseA, 
                                              # cmPat = cmPat))
                # }
                
                # RAND.EFF <- sapply(allow.doses, function(a){
                        # mean(apply(sim.alphas, 2, function(o){
                                # as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        # }))
                # })
                # RAND.EFF <- exp(RAND.EFF)/sum(exp(RAND.EFF))
                # nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)
                
                # ####a condition for untried higher dose level that is predicted to be efficacious####
                # if(nxtdose > max(rec_doseA) + 1 & !nxtdose %in% rec_doseA){
                        # nxtdose <- max(rec_doseA) + 1
                # }
        # }
        # #print('Stage 3');
        # ############################################################################
        # ###############################Stage 3######################################
        # ############################################################################
        # ###################when all the data are available##########################
        # ############################################################################
        # doseA <- nxtdose
        # rec_doseA <- c(rec_doseA, doseA)
        # altMat[doseA, k] <- altMat[doseA, k] + chSize
        # cmPat <- cmPat + chSize
        # PatData <- list()
        # n <- chSize
        # K <- MaxCycle
        
        # ##################################################
        # #######design matrix for toxicity#################
        # ##################################################
        # dose.pat <- rep(doseA, n)
        # xx <- expand.grid(dose.pat, 1:K)
        # X <- as.matrix(cbind(rep(1, n*K), xx))
        # X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
        # colnames(X) <- c("intercept", "dose", "cycle", "icycle")
        
        # ##################################################
        # #######design matrix for efficacy#################
        # ##################################################
        # Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
        # colnames(Xe) <- c("intercept", "dose", "squared dose")
        
        # #####simulate nTTP and Efficacy for the cohort########
        # outcome <- apply(X, 1, function(a){
                # tox.by.grade <- tox_matrix[a[2], a[3], ,]
                # nttp.indices <- apply(tox.by.grade, 1, function(o){
                        # sample(1:5, 1, prob = o)
                # })
                # if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                        # y.dlt <- 1
                # }else{
                        # y.dlt <- 0
                # }
                # y.nTTP <- nTTP.all[nttp.indices[1], 
                                   # nttp.indices[2],
                                   # nttp.indices[3]]
                # return(c(y.dlt = y.dlt,
                         # y = y.nTTP))
        # })
        # y <- outcome[2, ]
        # y.dlt <- outcome[1, ]
        
        # beta.mean <- eff.structure[Xe[1, 2]]
        # beta.sd <- eff.sd
        # beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
        # beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
        # e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
        # #################################################################
        # #################################################################
        
        # #####add to the master Dataset and extract from this dataset for final analysis###
        # #####PatData stores the dataset used to estimate the model########################
        # temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n*K), 1:n))
        # temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
        # masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
        
        # temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
        # masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
        
        # PatData$toxicity <- masterData$toxicity[, c(1, 2, 3, 4, 5, 6, 9)]
        # #####toxicity dropout########
        # #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
        # dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                # if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                        # return(rep(0, length(o)))
                # }else{
                        # o.temp <- rep(0, length(o))
                        # o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                        # return(o.temp)
                # }
        # }, simplify = FALSE)
        
        # dropout.dlt <- NULL
        
        # for(i in unique(PatData$toxicity$subID)){
                # dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
        # }
        
        # PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
        
        # #####efficacy dropout########
        # dropout.eff <- sapply(dlt.drop, function(a){
                # if(sum(a) == 0){
                        # return(0)
                # }else if(min(which(a == 1)) >= 4){
                        # return(0)
                # }else{
                        # return(1)
                # }
        # })
        
        # dropout.eff <- as.numeric(dropout.eff)
        # PatData$efficacy <- masterData$efficacy[dropout.eff == 0, ]
        # #################################################################
        # #################Model Estimation################################
        # #################################################################
        
        # nTTP <- PatData$toxicity[, 2]
        # X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "cycle")])
        # n.subj <- length(unique(PatData$toxicity[, "subID"]))
        # W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
        # W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                # which(a == unique(PatData$toxicity$subID))
        # })))] <- 1
        
        # EFF <- PatData$efficacy[, 1]
        # X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
        # keepeff.ind <- which(dropout.eff == 0)
        
        # model.file <- function()
        # {       
                # beta <- c(beta_other[1], beta_dose, beta_other[2])
                # for(k in 1:Nsub){
                        # u[k] ~ dnorm(0, tau_u)
                # }
                # for(i in 1:N1){
                        # y[i] ~ dnorm(mu[i], tau_e)
                        # mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                # }
                # for(j in 1:N2){
                        # e[j] ~ dnorm(mu_e[j], tau_f)
                        # mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                # }
                
                # beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                # beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                # alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                # gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                # tau_e ~ dunif(p1_tau_e, p2_tau_e)
                # tau_u ~ dunif(p1_tau_u, p2_tau_u)
                # tau_f ~ dunif(p1_tau_f, p2_tau_f)
        # }
        
        # mydata <- list(N1 = length(nTTP), N2 = length(EFF), Nsub = n.subj, 
                       # y = nTTP, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                       # W_y = W_y, X_e = X_e, p1_beta_other = c(0, 0), 
                       # p2_beta_other = diag(rep(0.001, 2)),
                       # p1_beta_dose = 0, p2_beta_dose = 1000,
                       # p1_alpha = c(0, 0, 0), p2_alpha = diag(rep(0.001, 3)), 
                       # p1_gamma0 = 0, p2_gamma0 = 0.001,
                       # p1_tau_e = 0, p2_tau_e = 1000, 
                       # p1_tau_u = 0, p2_tau_u = 1000, 
                       # p1_tau_f = 0, p2_tau_f = 1000)
        
        # path.model <- file.path(tempdir(), "model.file.txt")
        # R2WinBUGS::write.model(model.file, path.model)
        
        # inits.list <- list(list(beta_other = c(0.1, 0.1),
                                # beta_dose = 0.1,
                                # alpha = c(0.1, 0.1, 0.1),
                                # gamma0 = 0.1,
                                # tau_e = 0.1,
                                # tau_u = 0.1,
                                # tau_f = 0.1,
                                # .RNG.seed = sample(1:1e+06, size = 1), 
                                # .RNG.name = "base::Wichmann-Hill"))
        
        # jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                     # quiet = TRUE, inits = inits.list)
        
        # update(jagsobj, n.iter = n.iter, progress.bar = "none")
        # post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                       # "gamma0"), 
                                            # n.iter = n.iter, 
                                            # progress.bar = "none")
        
        # sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                     # post.samples$beta_other[2,,1]))
        
        # sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
        
        # ############redefine allowable doses#############
        # ####condition 1####
        # prob1.doses <- sapply(doses, function(a){
                # mean(apply(sim.betas, 2, function(o){
                        # as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd1)
                # }))
        # })
        
        # ####condition 2####
        # prob2.doses <- sapply(doses, function(a){
                # mean(apply(sim.betas, 2, function(o){
                        # as.numeric(mean(sapply(2:K, function(m){
                                # as.numeric(o[1] + o[2] * a + o[3] * m)
                        # })) <= thrd2)
                # }))
        # })
        
        # allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
        # if(length(allow.doses) == 0){
                # return(results = list(sc = sc, opt.dose = 0, 
                                      # altMat = altMat, 
                                      # masterData = masterData, 
                                      # stage1.allow = stage1.allow, 
                                      # dose_trace = rec_doseA, 
                                      # cmPat = cmPat))
        # }
        
        # effcy.doses <- sapply(allow.doses, function(a){
                # mean(apply(sim.alphas, 2, function(o){
                        # o[1] + o[2] * a + o[3] * a^2}))
        # })
        
        # proxy.eff.doses <- allow.doses[which(sapply(effcy.doses, function(a){
                # abs(a - max(effcy.doses)) <= proxy.thrd
        # }))]
        
        # ###recommend the lowest dose that is efficacious####
        # recom.doses <- min(proxy.eff.doses)
        # #print('output')
		
        # return(results = list(sc = sc, opt.dose = recom.doses, 
                              # altMat = altMat, masterData = masterData, 
                              # stage1.allow = stage1.allow, 
                              # effy.doses = proxy.eff.doses, 
                              # dose_trace = rec_doseA, 
                              # cmPat = cmPat))
# }

