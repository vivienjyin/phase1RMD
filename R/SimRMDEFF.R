###################################################################################
###########Simulate Data from the Matrix and Do Stage 1-3##########################
###################################################################################
###model: yij = beta0 + beta1xi + beta2tj + ui + epiij#############################
######### ei = alpha0 + alpha1xi + alpha2xi^2 + gamma0ui + epi2i###################
######### epi ~ N(0, sigma_e); ui ~ N(0, sigma_u); epi2 ~ N(0, sigma_f)############
###################################################################################
#rm(list = ls())
#library(rjags)
#library(R2WinBUGS)
#library(BiocParallel)
#library(arrayhelpers)


SimRMDEFF <- function(numTrials = 100, trialSize = 36, doses = 1:6, cycles = 1:6,  
                      eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9), 
                      eff.sd = 0.2, tox.target = 0.28, p1 = 0.2, p2 = 0.2, 
                      ps1 = 0.2,  StrDose = 1, chSize = 3, tox.matrix = NULL, 
                      proxy.thrd = 0.1, thrd1 = 0.28, thrd2 = 0.28, 
                      wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                    0, 0.5, 0.75, 1 , 1.5, 
                                    0, 0  , 0   , 0.5, 1  ), 
                                  byrow = T, ncol = 5), 
                      toxmax = 2.5, toxtype = NULL, intercept.alpha = NULL, 
                      coef.beta = NULL, cycle.gamma = NULL, 
                      param.ctrl = list()) {
        
        Trial.Simulation <- function(trialSize = 36, doses = 1:6, cycles = 1:6,  
                                     eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9), 
                                     eff.sd = 0.2, tox.target = 0.28, p1 = 0.2, p2 = 0.2, 
                                     ps1 = 0.2, seed, StrDose = 1, chSize = 3, tox.matrix = NULL, 
                                     proxy.thrd = 0.1, thrd1 = 0.28, thrd2 = 0.28, 
                                     wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                                   0, 0.5, 0.75, 1 , 1.5, 
                                                   0, 0  , 0   , 0.5, 1  ), 
                                                 byrow = T, ncol = 5), 
                                     toxmax = 2.5, toxtype = NULL, intercept.alpha = NULL, 
                                     coef.beta = NULL, cycle.gamma = NULL, 
                                     param.ctrl = list()){
                
                ctrl_param <- list(p1_beta_intercept = 0, 
                                   p1_beta_cycle = 0,
                                   p2_beta_intercept = 0.001,
                                   p2_beta_cycle = 0.001,
                                   p1_beta_dose = 0, 
                                   p2_beta_dose = 1000,
                                   p1_alpha = c(0, 0, 0), 
                                   p2_alpha = diag(rep(0.001, 3)), 
                                   p1_gamma0 = 0, 
                                   p2_gamma0 = 0.001)
                ctrl_param <- modifyList(ctrl_param, param.ctrl)
                
                listdo <- paste0("D", doses)
                if(length(eff.structure) != length(doses))
                        stop("Check if you have specified the efficacy mean structure for the correct number of doses.")
                if(is.null(wm))
                        stop("Make sure that you specified the clinical weight matrix for toxicities.")
                
                MaxCycle <- length(cycles)
                Numdose <- length(doses)
                flag.tox.matrix.null <- 0
                
                if(is.null(tox.matrix)){
                        flag.tox.matrix.null <- 1
                        tox.matrix <- array(NA, dim = c(Numdose, MaxCycle, nrow(wm), 5))
                        #if(length(toxtype) != 3)
                        #        stop("Right now we are only considering three toxicity types, but we will relax this constraint in the future work.")
                        if(length(intercept.alpha) != 4){
                                stop("Exactly four intercepts alpha are needed for grade 0--4 in simulation!")
                        } 
                        if(min(diff(intercept.alpha, lag = 1, differences = 1)) < 0){
                                stop("Intercepts alpha for simulation must be in a monotonic increasing order!")
                        }
                        #if(length(coef.beta) != 3){
                        #        stop("Exactly three betas need to be specified for three types of toxicities!")
                        #}
                }
                
                
                #############################################
                #########Determine Scenarios#################
                #############################################
                # w1 <- wm[1, ]
                # w2 <- wm[2, ]
                # w3 <- wm[3, ]
                
                ####Array of the normalized scores########
                # nTTP.array <- function(w1,w2,w3){
                #         nTTP <- array(NA,c(5,5,5))
                #         for (t1 in 1:5){
                #                 for (t2 in 1:5){
                #                         for (t3 in 1:5){
                #                                 nTTP[t1,t2,t3] <-  
                #                                         sqrt(w1[t1]^2 + w2[t2]^2 + w3[t3]^2)/toxmax
                #                         }
                #                 }
                #         }
                #         return(nTTP)
                # }
                # nTTP.all <- nTTP.array(w1,w2,w3)
                
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
                
                
                ####compute pdlt given probability array####
                pDLT <- function(proba){    
                        DLT_array <- array(NA, c(rep(5, nrow(wm))))
                        capacity <- 5^(nrow(wm))
                        for(tt in 1 : capacity) {
                                index <- vec2array(tt, dim = c(rep(5, nrow(wm))))
                                DLT_array[tt] <- as.numeric(max(wm[t(rbind(1 : nrow(wm), index))]) >= 1)
                        }
                        return(sum(proba * DLT_array))
                }
                
                ####compute mean nTTP given probability array######
                mnTTP <- function(proba){    
                        return(sum(proba * nTTP.all))
                }
                
                sc.mat <- matrix(NA, MaxCycle, Numdose)
                pdlt <- matrix(NA, MaxCycle, Numdose)
                
                for(i in 1:MaxCycle){   ##loop through cycle
                        for(j in 1:Numdose){ ##loop through dose
                                if(flag.tox.matrix.null == 0){
                                        nTTP.prob <- tox.matrix[j, i, ,]
                                }else{
                                        # initialize storing matrix
                                        cump <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cumulative probability
                                        celp <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cell probability
                                        
                                        for(k in 1:length(toxtype)){    
                                                # proportional odds model
                                                logitcp <- intercept.alpha + coef.beta[k] * j +  cycle.gamma * (i - 1)
                                                # cumulative probabilities
                                                cump[k, 1:4] <- exp(logitcp) / (1 + exp(logitcp)) 
                                                cump[k, 5] <- 1
                                                # cell probabilities
                                                celp[k, ] <- c(cump[k,1], diff(cump[k, ], lag = 1, differences = 1))
                                        }
                                        
                                        nTTP.prob <- celp
                                        tox.matrix[j, i, ,] <- celp
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
                                proba <- nTTP_prob(nTTP.prob)
                                sc.mat[i, j] <- round(mnTTP(proba), 4)
                                pdlt[i, j] <- round(pDLT(proba), 4)
                        }
                }
                
                mEFF <-  eff.structure
                sc <- rbind(sc.mat, mEFF, pdlt[1, ])
                colnames(sc) <- listdo
                rownames(sc) <- c("mnTTP.1st", "mnTTP.2nd", "mnTTP.3rd", 
                                  "mnTTP.4th", "mnTTP.5th", "mnTTP.6th",
                                  "mEFF", "pDLT")
                
                ####################################################################################
                set.seed(seed)
                
                k <- 1  ####Let us use k to denote the column of altMat####
                doseA <- StrDose ####Let us use doseA to denote the current dose assigned
                cmPat <- chSize ####current patient size
                masterData <- list()
                rec_doseA <- NULL
                stage1.allow <- NULL
                
                altMat <- matrix(0, Numdose, k)
                altMat[doseA, k] <- chSize
                
                ###############################################################
                ########################Stage 1################################
                ###############################################################
                
                while(cmPat <= trialSize/2){
                        rec_doseA <- c(rec_doseA, doseA)
                        PatData <- NULL
                        n <- chSize
                        K <- MaxCycle
                        
                        ##############################################################
                        #############design matrix for toxicity#######################
                        ##############################################################
                        dose.pat <- rep(doseA, n)
                        xx <- expand.grid(dose.pat, 1:K)
                        X <- as.matrix(cbind(rep(1, n * K), xx))
                        X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
                        colnames(X) <- c("intercept", "dose", "cycle", "icycle")
                        
                        ##############################################################
                        #######design matrix for efficacy#############################
                        ##############################################################
                        Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
                        colnames(Xe) <- c("intercept", "dose", "squared dose")
                        
                        #####simulate nTTP and Efficacy for the cohort########
                        outcome <- apply(X, 1, function(a){
                                tox.by.grade <- tox.matrix[a[2], a[3], ,]
                                nttp.indices <- apply(tox.by.grade, 1, function(o){
                                        sample(1:5, 1, prob = o)
                                })
                                if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                                        y.dlt <- 1
                                }else{
                                        y.dlt <- 0
                                }
                                #y.nTTP <- nTTP.all[nttp.indices[1], 
                                #                   nttp.indices[2],
                                #                   nttp.indices[3]]
                                y.nTTP <- nTTP.all[matrix(nttp.indices, nrow = 1)]
                                return(c(y.dlt = y.dlt,
                                         y = y.nTTP))
                        })
                        y <- outcome[2, ]
                        y.dlt <- outcome[1, ]
                        
                        beta.mean <- eff.structure[Xe[1, 2]]
                        beta.sd <- eff.sd
                        beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
                        beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
                        e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
                        
                        
                        #####construct master Dataset and extract from this dataset for each interim#########
                        #####PatData stores the dataset used to estimate the model at each interim###########
                        temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n * K), 1:n))
                        temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
                        masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
                        
                        temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
                        masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
                        
                        current_cohort <- cmPat/chSize
                        PatData.index <- data.frame(cohort = 1:current_cohort, cycles = current_cohort:1)
                        PatData.index <- within(PatData.index, cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles))
                        for(i in 1:nrow(PatData.index)){
                                PatData <- rbind(PatData, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                                      masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                              c(1, 2, 3, 4, 5, 6, 9)])
                        }
                        
                        #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt#######
                        dlt.drop <- sapply(by(PatData, PatData$subID, function(a){a[, 1]}), function(o){
                                if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                        return(rep(0, length(o)))
                                }else{
                                        o.temp <- rep(0, length(o))
                                        o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                        return(o.temp)
                                }
                        })
                        
                        dropout.dlt <- NULL
                        
                        for(i in unique(PatData$subID)){
                                dropout.dlt[PatData$subID == i] <- dlt.drop[[i]]
                        }
                        
                        PatData <- PatData[dropout.dlt == 0, ]
                        
                        ###################################################################################
                        ################Model Estimation given PatData#####################################
                        ################Stage 1 only considers the toxicity model##########################
                        ###################################################################################
                        nTTP <- PatData[, 2]
                        X_y <- as.matrix(PatData[, c("intercept", "dose", "icycle")])
                        n.subj <- length(unique(PatData[, "subID"]))
                        W_y <- matrix(0, nrow(PatData), n.subj)
                        W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$subID, function(a){
                                which(a == unique(PatData$subID))
                        })))] <- 1
						beta_other <- 1; beta_dose <- 1; N1 <- 1; inprod <- function(){}; u <- 1;
						N2 <- 1; Dose <- 1; Post_nTTP <- 1; Nsub <- 1; gamma0 <- 1;
						Post_EFF <- 1;                        
                        model.file <- function()
                        {       
                                beta <- c(beta_other[1], beta_dose, beta_other[2])
                                for(i in 1:N1){
                                        y[i] ~ dnorm(mu[i], tau_e)
                                        mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                                }
                                for(j in 1:N2){
                                        u[j] ~ dnorm(0, tau_u)
                                }
                                beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                                beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                                tau_e ~ dunif(p1_tau_e, p2_tau_e)
                                tau_u ~ dunif(p1_tau_u, p2_tau_u)
                        }
                        
                        mydata <- list(N1 = length(nTTP), N2 = n.subj, y = nTTP, X_y = X_y, 
                                       W_y = W_y, p1_beta_other = c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle), 
                                       p2_beta_other = diag(c(ctrl_param$p2_beta_intercept, ctrl_param$p2_beta_cycle)),
                                       p1_beta_dose = ctrl_param$p1_beta_dose, p2_beta_dose = ctrl_param$p2_beta_dose,
                                       p1_tau_e = 0, p2_tau_e = 1000, 
                                       p1_tau_u = 0, p2_tau_u = 1000)
                        
                        path.model <- file.path(tempdir(), "model.file.txt")
                        #R2WinBUGS::write.model(model.file, path.model)
                        
                        inits.list <- list(list(beta_other = c(0.1, 0.1),
                                                beta_dose = 0.1,
                                                tau_e = 0.1,
                                                tau_u = 0.1,
                                                .RNG.seed = sample(1:1e+06, size = 1), 
                                                .RNG.name = "base::Wichmann-Hill"))
                        
                        jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                                     quiet = TRUE, inits = inits.list)
                        
                        update(jagsobj, n.iter = 5000, progress.bar = "none")
                        post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other"), 
                                                            n.iter = 5000, 
                                                            progress.bar = "none")
                        
                        if(cmPat == trialSize/2){
                                ######################################################################
                                #############define allowable doses###################################
                                ######################################################################
                                sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                             post.samples$beta_other[2,,1]))
                                ####condition 1####
                                prob1.doses <- sapply(doses, function(a){
                                        mean(apply(sim.betas, 2, function(o){
                                                as.numeric(o[1] + o[2] * a <= thrd1)
                                        }))
                                })
                                
                                ####condition 2####
                                prob2.doses <- sapply(doses, function(a){
                                        mean(apply(sim.betas, 2, function(o){
                                                as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                                        }))
                                })
                                
                                allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                                if(length(allow.doses) == 0){
                                        stage1.allow <- 0
                                        return(results = list(sc = sc, opt.dose = 0, 
                                                              altMat = altMat, 
                                                              masterData = masterData, 
                                                              stage1.allow = stage1.allow, 
                                                              dose_trace = rec_doseA, 
                                                              cmPat = cmPat))
                                }
                                stage1.allow <- allow.doses
                                
                        }else{
                                sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1]))
                                loss.doses <- sapply(doses, function(a){
                                        mean(apply(sim.betas, 2, function(o){
                                                abs(o[1] + o[2] * a - tox.target)
                                        }))
                                })
                                nxtdose <- doses[which.min(loss.doses)]
                                
                                if(as.numeric(nxtdose) > (doseA + 1)){
                                        doseA <- doseA + 1;
                                }else{
                                        doseA <- as.numeric(nxtdose);
                                }
                                
                                if (cmPat == chSize){
                                        dlt <- with(masterData$toxicity, y.dlt[cycle == 1 & V7 == cmPat/chSize])
                                        if (sum(dlt) == 0){
                                                doseA <- 2
                                                dose_flag <- 0
                                        }else{
                                                doseA <- 1
                                                dose_flag <- 1			
                                        }
                                }
                                if ((cmPat >= chSize*2) & (dose_flag == 1)){
                                        dlt <- with(masterData$toxicity, y.dlt[cycle == 1 & V7 == cmPat/chSize])
                                        if (sum(dlt) == 0){
                                                doseA <- 2
                                                dose_flag <- 0
                                        }else{
                                                doseA <- 1
                                                dose_flag <- 1			
                                        }
                                }
                                if(cmPat >= chSize*3){
                                        sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                                     post.samples$beta_other[2,,1]))
                                        ####condition 1####
                                        prob1.doses <- sapply(doses, function(a){
                                                mean(apply(sim.betas, 2, function(o){
                                                        as.numeric(o[1] + o[2] * a <= thrd1)
                                                }))
                                        })
                                        
                                        ####condition 2####
                                        prob2.doses <- sapply(doses, function(a){
                                                mean(apply(sim.betas, 2, function(o){
                                                        as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                                                }))
                                        })
                                        
                                        allow.doses <- which(prob1.doses >= ps1 & prob2.doses >= ps1)
                                        if(length(allow.doses) == 0){
                                                stage1.allow <- 0
                                                return(results = list(sc = sc, opt.dose = 0, 
                                                                      altMat = altMat, 
                                                                      masterData = masterData, 
                                                                      stage1.allow = stage1.allow, 
                                                                      dose_trace = rec_doseA, 
                                                                      cmPat = cmPat))
                                        }
                                }
                                altMat[doseA, k] <- altMat[doseA, k] + chSize
                        }
                        cmPat <- cmPat + chSize # cmPat for the next cohort (loop)
                }
                
                ###############################################################
                ########################Stage 2################################
                ###############################################################
                #######let us wait until the efficacy data for stage 1 is available######
                #########################################################################
                cmPat <- cmPat - chSize
                
                
                PatData <- list() 
                PatData.index <- data.frame(cohort = current_cohort:1, cycles = 3:(3 + current_cohort - 1))
                PatData.index <- within(PatData.index, cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles))
                for(i in 1:nrow(PatData.index)){
                        PatData$toxicity <- rbind(PatData$toxicity, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                                                masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                                        c(1, 2, 3, 4, 5, 6, 9)])
                }
                
                #####change the ordering of the subjects--being consistent######
                temp.toxicity <- NULL
                for(i in 1:current_cohort){
                        temp.toxicity <- rbind(temp.toxicity, PatData$toxicity[grepl(paste0("cohort",i,"subject"), PatData$toxicity$subID), ])
                }
                
                PatData$toxicity <- temp.toxicity
                rm(temp.toxicity)
                
                #####toxicity dropout########
                #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
                dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                        if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                return(rep(0, length(o)))
                        }else{
                                o.temp <- rep(0, length(o))
                                o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                return(o.temp)
                        }
                })
                
                dropout.dlt <- NULL
                
                for(i in unique(PatData$toxicity$subID)){
                        dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
                }
                
                PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
                
                #####efficacy dropout########
                dropout.eff <- sapply(dlt.drop, function(a){
                        if(sum(a) == 0){
                                return(0)
                        }else if(min(which(a == 1)) >= 4){
                                return(0)
                        }else{
                                return(1)
                        }
                })
                
                dropout.eff <- as.numeric(dropout.eff)
                PatData$efficacy <- masterData$efficacy[dropout.eff == 0, ]
                
                
                ############Model Estimation to randomize among allowable doses##########
                nTTP <- PatData$toxicity[, 2]
                X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
                n.subj <- length(unique(PatData$toxicity[, "subID"]))
                W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                        which(a == unique(PatData$toxicity$subID))
                })))] <- 1
                
                EFF <- PatData$efficacy[, 1]
                X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
                keepeff.ind <- which(dropout.eff == 0)
                
                model.file <- function()
                {       
                        beta <- c(beta_other[1], beta_dose, beta_other[2])
                        for(i in 1:N1){
                                y[i] ~ dnorm(mu[i], tau_e)
                                mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                        }
                        for(k in 1:Nsub){
                                u[k] ~ dnorm(0, tau_u)
                        }
                        for(j in 1:N2){
                                e[j] ~ dnorm(mu_e[j], tau_f)
                                mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                        }
                        beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                        beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                        alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                        gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                        tau_e ~ dunif(p1_tau_e, p2_tau_e)
                        tau_u ~ dunif(p1_tau_u, p2_tau_u)
                        tau_f ~ dunif(p1_tau_f, p2_tau_f)
                }
                
                mydata <- list(N1 = length(nTTP), N2 = length(EFF), y = nTTP, 
                               Nsub = n.subj, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                               W_y = W_y, X_e = X_e, p1_beta_other = c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle), 
                               p2_beta_other = diag(c(ctrl_param$p2_beta_intercept, ctrl_param$p2_beta_cycle)),
                               p1_beta_dose = ctrl_param$p1_beta_dose, p2_beta_dose = ctrl_param$p2_beta_dose,
                               p1_alpha = ctrl_param$p1_alpha, p2_alpha = ctrl_param$p2_alpha, 
                               p1_gamma0 = ctrl_param$p1_gamma0, p2_gamma0 = ctrl_param$p2_gamma0,
                               p1_tau_e = 0, p2_tau_e = 1000, 
                               p1_tau_u = 0, p2_tau_u = 1000, 
                               p1_tau_f = 0, p2_tau_f = 1000)
                
                path.model <- file.path(tempdir(), "model.file.txt")
                #R2WinBUGS::write.model(model.file, path.model)
                
                inits.list <- list(list(beta_other = c(0.1, 0.1),
                                        beta_dose = 0.1,
                                        alpha = c(0.1, 0.1, 0.1),
                                        gamma0 = 0.1,
                                        tau_e = 0.1,
                                        tau_u = 0.1,
                                        tau_f = 0.1,
                                        .RNG.seed = sample(1:1e+06, size = 1), 
                                        .RNG.name = "base::Wichmann-Hill"))
                
                jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                             quiet = TRUE, inits = inits.list)
                
                update(jagsobj, n.iter = 5000, progress.bar = "none")
                post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                               "gamma0"), 
                                                    n.iter = 5000, 
                                                    progress.bar = "none")
                
                ######################################################################
                #############update allowable doses###################################
                ######################################################################
                sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                             post.samples$beta_other[2,,1]))
                ####condition 1####
                prob1.doses <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a <= thrd1)
                        }))
                })
                
                ####condition 2####
                prob2.doses <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                        }))
                })
                
                allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                if(length(allow.doses) == 0){
                        return(results = list(sc = sc, opt.dose = 0, 
                                              altMat = altMat, 
                                              masterData = masterData, 
                                              stage1.allow = stage1.allow, 
                                              dose_trace = rec_doseA, 
                                              cmPat = cmPat))
                }
                
                sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                RAND.EFF <- sapply(allow.doses, function(a){
                        mean(apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        }))
                })
                RAND.EFF <- exp(RAND.EFF)/sum(exp(RAND.EFF))
                nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)
                
                ####a condition for untried higher dose level that is randomized#### rec_doseA[length(rec_doseA)]
                if(nxtdose > max(rec_doseA) + 1 & !nxtdose %in% rec_doseA){
                        nxtdose <- max(rec_doseA) + 1
                }
                
                
                #################################################################
                ########generate data for the new enrolled cohort in Stage 2#####
                #################################################################
                while(cmPat < trialSize - chSize){
                        doseA <- nxtdose
                        rec_doseA <- c(rec_doseA, doseA)
                        altMat[doseA, k] <- altMat[doseA, k] + chSize
                        cmPat <- cmPat + chSize
                        PatData <- list()
                        n <- chSize
                        K <- MaxCycle
                        
                        #######################################################
                        #################design matrix for toxicity############
                        #######################################################
                        dose.pat <- rep(doseA, n)
                        xx <- expand.grid(dose.pat, 1:K)
                        X <- as.matrix(cbind(rep(1, n*K), xx))
                        X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
                        colnames(X) <- c("intercept", "dose", "cycle", "icycle")
                        
                        #########################################################
                        #######design matrix for efficacy########################
                        #########################################################
                        Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
                        colnames(Xe) <- c("intercept", "dose", "squared dose")
                        
                        #####simulate nTTP and Efficacy for the cohort########
                        outcome <- apply(X, 1, function(a){
                                tox.by.grade <- tox.matrix[a[2], a[3], ,]
                                nttp.indices <- apply(tox.by.grade, 1, function(o){
                                        sample(1:5, 1, prob = o)
                                })
                                if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                                        y.dlt <- 1
                                }else{
                                        y.dlt <- 0
                                }
                                # y.nTTP <- nTTP.all[nttp.indices[1], 
                                #                    nttp.indices[2],
                                #                    nttp.indices[3]]
                                y.nTTP <- nTTP.all[matrix(nttp.indices, nrow = 1)]
                                return(c(y.dlt = y.dlt,
                                         y = y.nTTP))
                        })
                        y <- outcome[2, ]
                        y.dlt <- outcome[1, ]
                        
                        beta.mean <- eff.structure[Xe[1, 2]]
                        beta.sd <- eff.sd
                        beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
                        beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
                        e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
                        #################################################################
                        #################################################################
                        
                        #####add to the master Dataset and extract from this dataset for each interim#########
                        #####PatData stores the dataset used to estimate the model at each interim############
                        temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n*K), 1:n))
                        temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
                        masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
                        
                        temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
                        masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
                        
                        current_cohort <- cmPat/chSize
                        PatData.index <- data.frame(cohort = 1:current_cohort, cycles = current_cohort:1)
                        cycles.adj <- c(rep(2, trialSize/(2 * chSize)), rep(0, current_cohort - trialSize/(2 * chSize)))
                        PatData.index <- within(PatData.index, {cycles <- cycles + cycles.adj
                        cycles <- ifelse(cycles >= MaxCycle, MaxCycle, cycles)})
                        for(i in 1:nrow(PatData.index)){
                                PatData$toxicity <- rbind(PatData$toxicity, masterData$toxicity[masterData$toxicity[, 7] == PatData.index[i, "cohort"] & 
                                                                                                        masterData$toxicity[, 5] <= PatData.index[i, "cycles"], 
                                                                                                c(1, 2, 3, 4, 5, 6, 9)])
                        }
                        
                        #####toxicity dropout########
                        #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
                        dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                                if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                        return(rep(0, length(o)))
                                }else{
                                        o.temp <- rep(0, length(o))
                                        o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                        return(o.temp)
                                }
                        })
                        
                        dropout.dlt <- NULL
                        for(i in unique(PatData$toxicity$subID)){
                                dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
                        }
                        
                        PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
                        cohort.eff.index <- with(PatData.index, which(cycles >= 3))
                        
                        #####efficacy dropout########
                        dropout.eff <- sapply(dlt.drop, function(a){
                                if(sum(a) == 0){
                                        return(0)
                                }else if(min(which(a == 1)) >= 4){
                                        return(0)
                                }else{
                                        return(1)
                                }
                        })
                        
                        dropout.eff <- as.numeric(dropout.eff)
                        
                        PatData$efficacy <- subset(masterData$efficacy, masterData$efficacy[, 5] %in% cohort.eff.index & dropout.eff == 0)
                        
                        
                        
                        #############################################################################################################################
                        ######################Model Estimation to update toxicity information and randomization probs################################
                        #############################################################################################################################
                        nTTP <- PatData$toxicity[, 2]
                        X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
                        n.subj <- length(unique(PatData$toxicity[, "subID"]))
                        W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                        W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                                which(a == unique(PatData$toxicity$subID))
                        })))] <- 1
                        
                        EFF <- PatData$efficacy[, 1]
                        X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
                        keepeff.ind <- which(dropout.eff == 0 & masterData$efficacy[, 5] %in% cohort.eff.index)
                        
                        model.file <- function()
                        {       
                                beta <- c(beta_other[1], beta_dose, beta_other[2])
                                for(k in 1:Nsub){
                                        u[k] ~ dnorm(0, tau_u)
                                }
                                for(i in 1:N1){
                                        y[i] ~ dnorm(mu[i], tau_e)
                                        mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                                }
                                for(j in 1:N2){
                                        e[j] ~ dnorm(mu_e[j], tau_f)
                                        mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                                }
                                
                                beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                                beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                                alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                                gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                                tau_e ~ dunif(p1_tau_e, p2_tau_e)
                                tau_u ~ dunif(p1_tau_u, p2_tau_u)
                                tau_f ~ dunif(p1_tau_f, p2_tau_f)
                        }
                        
                        mydata <- list(N1 = length(nTTP), N2 = length(EFF), y = nTTP, 
                                       Nsub = n.subj, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                                       W_y = W_y, X_e = X_e, p1_beta_other = c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle), 
                                       p2_beta_other = diag(c(ctrl_param$p2_beta_intercept, ctrl_param$p2_beta_cycle)),
                                       p1_beta_dose = ctrl_param$p1_beta_dose, p2_beta_dose = ctrl_param$p2_beta_dose,
                                       p1_alpha = ctrl_param$p1_alpha, p2_alpha = ctrl_param$p2_alpha, 
                                       p1_gamma0 = ctrl_param$p1_gamma0, p2_gamma0 = ctrl_param$p2_gamma0,
                                       p1_tau_e = 0, p2_tau_e = 1000, 
                                       p1_tau_u = 0, p2_tau_u = 1000, 
                                       p1_tau_f = 0, p2_tau_f = 1000)
                        
                        path.model <- file.path(tempdir(), "model.file.txt")
                        #R2WinBUGS::write.model(model.file, path.model)
                        
                        inits.list <- list(list(beta_other = c(0.1, 0.1),
                                                beta_dose = 0.1,
                                                alpha = c(0.1, 0.1, 0.1),
                                                gamma0 = 0.1,
                                                tau_e = 0.1,
                                                tau_u = 0.1,
                                                tau_f = 0.1,
                                                .RNG.seed = sample(1:1e+06, size = 1), 
                                                .RNG.name = "base::Wichmann-Hill"))
                        
                        jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                                     quiet = TRUE, inits = inits.list)
                        
                        update(jagsobj, n.iter = 5000, progress.bar = "none")
                        post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                                       "gamma0"), 
                                                            n.iter = 5000, 
                                                            progress.bar = "none")
                        
                        sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                     post.samples$beta_other[2,,1]))
                        
                        sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                        
                        ############redefine allowable doses#############
                        ####condition 1####
                        prob1.doses <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a <= thrd1)
                                }))
                        })
                        
                        ####condition 2####
                        prob2.doses <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd2)
                                }))
                        })
                        
                        allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                        if(length(allow.doses) == 0){
                                return(results = list(sc = sc, opt.dose = 0, 
                                                      altMat = altMat, 
                                                      masterData = masterData, 
                                                      stage1.allow = stage1.allow, 
                                                      dose_trace = rec_doseA, 
                                                      cmPat = cmPat))
                        }
                        
                        RAND.EFF <- sapply(allow.doses, function(a){
                                mean(apply(sim.alphas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3] * a^2)
                                }))
                        })
                        RAND.EFF <- exp(RAND.EFF)/sum(exp(RAND.EFF))
                        nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)
                        
                        ####a condition for untried higher dose level that is predicted to be efficacious####
                        if(nxtdose > max(rec_doseA) + 1 & !nxtdose %in% rec_doseA){
                                nxtdose <- max(rec_doseA) + 1
                        }
                }
                
                ############################################################################
                ###############################Stage 3######################################
                ############################################################################
                ###################when all the data are available##########################
                ############################################################################
                doseA <- nxtdose
                rec_doseA <- c(rec_doseA, doseA)
                altMat[doseA, k] <- altMat[doseA, k] + chSize
                cmPat <- cmPat + chSize
                PatData <- list()
                n <- chSize
                K <- MaxCycle
                
                ##################################################
                #######design matrix for toxicity#################
                ##################################################
                dose.pat <- rep(doseA, n)
                xx <- expand.grid(dose.pat, 1:K)
                X <- as.matrix(cbind(rep(1, n*K), xx))
                X <- cbind(X, ifelse(X[, 3] == 1, 0, 1))
                colnames(X) <- c("intercept", "dose", "cycle", "icycle")
                
                ##################################################
                #######design matrix for efficacy#################
                ##################################################
                Xe <- cbind(1, rep(doseA, n), rep(doseA^2, n))
                colnames(Xe) <- c("intercept", "dose", "squared dose")
                
                #####simulate nTTP and Efficacy for the cohort########
                outcome <- apply(X, 1, function(a){
                        tox.by.grade <- tox.matrix[a[2], a[3], ,]
                        nttp.indices <- apply(tox.by.grade, 1, function(o){
                                sample(1:5, 1, prob = o)
                        })
                        if(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1){
                                y.dlt <- 1
                        }else{
                                y.dlt <- 0
                        }
                        # y.nTTP <- nTTP.all[nttp.indices[1], 
                        #                    nttp.indices[2],
                        #                    nttp.indices[3]]
                        y.nTTP <- nTTP.all[matrix(nttp.indices, nrow = 1)]
                        return(c(y.dlt = y.dlt,
                                 y = y.nTTP))
                })
                y <- outcome[2, ]
                y.dlt <- outcome[1, ]
                
                beta.mean <- eff.structure[Xe[1, 2]]
                beta.sd <- eff.sd
                beta.shape1 <- ((1 - beta.mean)/beta.sd^2 - 1/beta.mean) * beta.mean^2
                beta.shape2 <- beta.shape1 * (1/beta.mean - 1)
                e <- matrix(rbeta(chSize, beta.shape1, beta.shape2), ncol = 1)
                #################################################################
                #################################################################
                
                #####add to the master Dataset and extract from this dataset for final analysis###
                #####PatData stores the dataset used to estimate the model########################
                temp.mtdata <- data.frame(cbind(y.dlt, y, X, rep(cmPat/chSize, n*K), 1:n))
                temp.mtdata$subID <- paste0("cohort", temp.mtdata[, 7], "subject", temp.mtdata[, 8])
                masterData$toxicity <- data.frame(rbind(masterData$toxicity, temp.mtdata))
                
                temp.eff.mtdata <- data.frame(cbind(e, Xe, cmPat/chSize))
                masterData$efficacy <- data.frame(rbind(masterData$efficacy, temp.eff.mtdata))
                
                PatData$toxicity <- masterData$toxicity[, c(1, 2, 3, 4, 5, 6, 9)]
                #####toxicity dropout########
                #####dropout.dlt is 0 means to keep the observation; 1 means to drop it due to dlt
                dlt.drop <- sapply(by(PatData$toxicity, PatData$toxicity$subID, function(a){a[, 1]}), function(o){
                        if((length(o) > 1 & all(o[1: (length(o) - 1)] == 0)) | length(o) == 1){
                                return(rep(0, length(o)))
                        }else{
                                o.temp <- rep(0, length(o))
                                o.temp[(min(which(o == 1)) + 1) : length(o.temp)] <- 1
                                return(o.temp)
                        }
                }, simplify = FALSE)
                
                dropout.dlt <- NULL
                
                for(i in unique(PatData$toxicity$subID)){
                        dropout.dlt[PatData$toxicity$subID == i] <- dlt.drop[[i]]
                }
                
                PatData$toxicity <- PatData$toxicity[dropout.dlt == 0, ]
                
                #####efficacy dropout########
                dropout.eff <- sapply(dlt.drop, function(a){
                        if(sum(a) == 0){
                                return(0)
                        }else if(min(which(a == 1)) >= 4){
                                return(0)
                        }else{
                                return(1)
                        }
                })
                
                dropout.eff <- as.numeric(dropout.eff)
                PatData$efficacy <- masterData$efficacy[dropout.eff == 0, ]
                #################################################################
                #################Model Estimation################################
                #################################################################
                
                nTTP <- PatData$toxicity[, 2]
                X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "cycle")])
                n.subj <- length(unique(PatData$toxicity[, "subID"]))
                W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                        which(a == unique(PatData$toxicity$subID))
                })))] <- 1
                
                EFF <- PatData$efficacy[, 1]
                X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "squared.dose")])
                keepeff.ind <- which(dropout.eff == 0)
                
                model.file <- function()
                {       
                        beta <- c(beta_other[1], beta_dose, beta_other[2])
                        for(k in 1:Nsub){
                                u[k] ~ dnorm(0, tau_u)
                        }
                        for(i in 1:N1){
                                y[i] ~ dnorm(mu[i], tau_e)
                                mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], u)
                        }
                        for(j in 1:N2){
                                e[j] ~ dnorm(mu_e[j], tau_f)
                                mu_e[j] <- inprod(X_e[j, ], alpha) + gamma0 * u[keepeff.ind[j]]
                        }
                        
                        beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
                        beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
                        alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
                        gamma0 ~ dnorm(p1_gamma0, p2_gamma0)
                        tau_e ~ dunif(p1_tau_e, p2_tau_e)
                        tau_u ~ dunif(p1_tau_u, p2_tau_u)
                        tau_f ~ dunif(p1_tau_f, p2_tau_f)
                }
                
                mydata <- list(N1 = length(nTTP), N2 = length(EFF), y = nTTP, 
                               Nsub = n.subj, X_y = X_y, e = EFF, keepeff.ind = keepeff.ind, 
                               W_y = W_y, X_e = X_e, p1_beta_other = c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle), 
                               p2_beta_other = diag(c(ctrl_param$p2_beta_intercept, ctrl_param$p2_beta_cycle)),
                               p1_beta_dose = ctrl_param$p1_beta_dose, p2_beta_dose = ctrl_param$p2_beta_dose,
                               p1_alpha = ctrl_param$p1_alpha, p2_alpha = ctrl_param$p2_alpha, 
                               p1_gamma0 = ctrl_param$p1_gamma0, p2_gamma0 = ctrl_param$p2_gamma0,
                               p1_tau_e = 0, p2_tau_e = 1000, 
                               p1_tau_u = 0, p2_tau_u = 1000, 
                               p1_tau_f = 0, p2_tau_f = 1000)
                
                path.model <- file.path(tempdir(), "model.file.txt")
                #R2WinBUGS::write.model(model.file, path.model)
                
                inits.list <- list(list(beta_other = c(0.1, 0.1),
                                        beta_dose = 0.1,
                                        alpha = c(0.1, 0.1, 0.1),
                                        gamma0 = 0.1,
                                        tau_e = 0.1,
                                        tau_u = 0.1,
                                        tau_f = 0.1,
                                        .RNG.seed = sample(1:1e+06, size = 1), 
                                        .RNG.name = "base::Wichmann-Hill"))
                
                jagsobj <- rjags::jags.model(path.model, data = mydata, n.chains = 1, 
                                             quiet = TRUE, inits = inits.list)
                
                update(jagsobj, n.iter = 5000, progress.bar = "none")
                post.samples <- rjags::jags.samples(jagsobj, c("beta_dose", "beta_other", "alpha",
                                                               "gamma0"), 
                                                    n.iter = 5000, 
                                                    progress.bar = "none")
                
                sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                             post.samples$beta_other[2,,1]))
                
                sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                
                ############redefine allowable doses#############
                ####condition 1####
                prob1.doses <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * 1 <= thrd1)
                        }))
                })
                
                ####condition 2####
                prob2.doses <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(mean(sapply(2:K, function(m){
                                        as.numeric(o[1] + o[2] * a + o[3] * m)
                                })) <= thrd2)
                        }))
                })
                
                allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                if(length(allow.doses) == 0){
                        return(results = list(sc = sc, opt.dose = 0, 
                                              altMat = altMat, 
                                              masterData = masterData, 
                                              stage1.allow = stage1.allow, 
                                              dose_trace = rec_doseA, 
                                              cmPat = cmPat))
                }
                
                effcy.doses <- sapply(allow.doses, function(a){
                        mean(apply(sim.alphas, 2, function(o){
                                o[1] + o[2] * a + o[3] * a^2}))
                })
                
                proxy.eff.doses <- allow.doses[which(sapply(effcy.doses, function(a){
                        abs(a - max(effcy.doses)) <= proxy.thrd
                }))]
                
                ###recommend the lowest dose that is efficacious####
                recom.doses <- min(proxy.eff.doses)
                
                return(results = list(sc = sc, opt.dose = recom.doses, 
                                      altMat = altMat, masterData = masterData, 
                                      stage1.allow = stage1.allow, 
                                      effy.doses = proxy.eff.doses, 
                                      dose_trace = rec_doseA, 
                                      cmPat = cmPat))
        }
        
        
        results.summary <- function(result){
                ########recommendation percentage###########
                ind <- NULL
                opt.doses <- sapply(result, function(a){
                        a$opt.dose
                })
                if(class(opt.doses) %in% c("integer", "numeric")){
                        recom.perc <- sapply(c(0, doses), function(a){sum(opt.doses == a)/length(opt.doses)})
                }else{
                        ind <- which(sapply(opt.doses, function(a){
                                is.null(a)
                        }))
                        opt.temp <- as.vector(do.call(rbind, opt.doses[setdiff(1:length(result), ind)]))
                        recom.perc <- sapply(c(0, doses), function(a){sum(opt.temp == a)/length(opt.temp)})
                }
                
                ########allocation percentage##################
                alloc <- sapply(result, function(a){
                        a$altMat
                })
                if(class(alloc) == "list"){
                        alloc.data <- do.call(cbind, alloc[setdiff(1:length(result), ind)])
                        alloc.perc <- apply(alloc.data, 1, function(a){
                                sum(a)/sum(alloc.data)
                        })
                }else{
                        alloc.perc <- apply(alloc, 1, function(a){
                                sum(a)/sum(alloc)
                        })
                }
                
                ###########stage 1 allowable doses region percentage##############
                stage1.allowed <- sapply(result, function(a){
                        a$stage1.allow
                })
                stage1.allowed <- stage1.allowed[setdiff(1:length(result), ind)]
                stage1.allow.perc <- sapply(c(0, doses), function(a){
                        mean(sapply(stage1.allowed, function(o){
                                as.numeric(a %in% o)
                        }))
                })
                
                ###########efficacious doses perc##################################
                effy.doses <- sapply(result, function(a){
                        a$effy.doses
                })
                effy.doses <- effy.doses[setdiff(1:length(result), ind)]
                effy.perc <- sapply(doses, function(a){
                        mean(sapply(effy.doses, function(o){
                                as.numeric(a %in% o)
                        }))
                })
                effy.null <- sum(sapply(effy.doses, function(a){
                        as.numeric(is.null(a))
                }))/length(effy.doses)
                
                ########scenario###########
                i <- 1
                while(is.null(result[[i]]$sc)){
                        i <- i + 1
                }
                sc <- result[[i]]$sc
                
                op.table <- rbind(round(recom.perc[2:(length(doses) + 1)], 3), round(alloc.perc, 3), 
                                  round(stage1.allow.perc[2:(length(doses) + 1)], 3), round(effy.perc, 3))
                op.table <- cbind(op.table, NA)
                op.table[c(1, 3, 4), (length(doses) + 1)] <- c(round(recom.perc[1], 3), 
                                                               round(stage1.allow.perc[1], 3), 
                                                               round(effy.null, 3))
                
                row.names(op.table) <- c("Recommendation (%)", 
                                         "Allocation (%)", 
                                         "Allowable Doses_Stage 1 (%)", 
                                         "Efficacious Doses (%)")
                colnames(op.table) <- c(paste0("D", doses), "None")
                
                
                # result.table <- rbind(sc, round(recom.perc[2:7], 4), round(alloc.perc, 4), 
                #                       round(stage1.allow.perc[2:7], 4), round(effy.perc, 4))
                # result.table <- cbind(result.table, "None")
                # row.names(result.table) <- c("mnTTP.1st", "mnTTP.2nd", "mnTTP.3rd", 
                #                              "mnTTP.4th", "mnTTP.5th", "mnTTP.6th", "mEFF", "pDLT", 
                #                              "recom.perc", "alloc.perc", "stage1allow.perc", 
                #                              "effy.region.perc")
                # result.table[c(9, 11, 12), 7] <- c(round(recom.perc[1], 4), 
                #                                    round(stage1.allow.perc[1], 4), 
                #                                    round(effy.null, 4))
                # result.table[result.table == "None"] <- ""
                # 
                # #####average patients enrolled #####
                # no.pat <- sapply(result, function(a){
                #         if(is.null(a$effy.nonmissing)){
                #                 c(a$cmPat, 0)
                #         }else{
                #                 c(a$cmPat, nrow(a$effy.nonmissing))
                #         }
                #         
                # })
                # if(class(no.pat) %in% c("integer", "numeric", "matrix")){
                #         avg.no.pat <- mean(no.pat[1, ])
                #         effy.nonm.perc <- mean(no.pat[2, ]/no.pat[1, ])
                # }else{
                #         no.pat.temp <- do.call(rbind, no.pat[setdiff(1:length(result), ind)])
                #         avg.no.pat <- mean(no.pat.temp[, 1])
                #         effy.nonm.perc <- mean(no.pat.temp[, 2]/no.pat.temp[, 1])
                # }
                # 
                
                #return(list(table.output = result.table, cmPat = avg.no.pat, effy.nonm.perc = effy.nonm.perc))
                return(list(op.table = op.table, sc = sc))        
                
        }
        set.seed(numTrials + 3)
        list_simul <- list()
        seed_rand <- sample(1 : 1000000000, numTrials)
        for(u in seq_along(seed_rand)) {
                list_simul[[u]] <- Trial.Simulation(trialSize = trialSize, doses = doses, cycles = cycles,  
                                                    eff.structure = eff.structure, 
                                                    eff.sd = eff.sd, tox.target = tox.target, p1 = p1, p2 = p2, 
                                                    ps1 = ps1, seed = seed_rand[u], StrDose = StrDose, chSize = chSize, tox.matrix = tox.matrix, 
                                                    proxy.thrd = proxy.thrd, thrd1 = thrd1, thrd2 = thrd2, 
                                                    wm = wm, toxmax = toxmax, toxtype = toxtype, 
                                                    intercept.alpha = intercept.alpha, 
                                                    coef.beta = coef.beta, cycle.gamma = cycle.gamma, 
                                                    param.ctrl = param.ctrl)
        }
        
        cat(sprintf("Operating characteristics based on %d simulations\n\n\nSample Size: %d\n\n\n", numTrials, trialSize))
        
        return(results.summary(list_simul))
        
}
