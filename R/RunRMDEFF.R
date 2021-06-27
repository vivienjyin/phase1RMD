#rm(list = ls())
#library(rjags)
#library(R2WinBUGS)
#library(ggplot2)

RunRMDEFF <- function(efficacy.dat = NULL, toxicity.dat, trialSize = 36,
                      seed = 624, chSize = 3, MaxCycle = 6, doses = 1:6, 
                      tox.target = 0.28, p1 = 0.2, p2 = 0.2, ps1 = 0.2, 
                      thrd1 = 0.28, thrd2 = 0.28, proxy.thrd = 0.1, 
                      dose_flag = 0, param.ctrl = list()){
        set.seed(seed)
        
        cat("Model : RMD with longitudinal toxicity\n\n\n")
        cat("Doses(skeleton):\n", paste0("\t", doses), "\n\n\n")
        
        subID <- toxicity.dat[, "subID"]
        cmPat <- max(as.numeric(do.call(rbind, strsplit(subID, "cohort|subject"))[, 2])) * 3
        
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
        PatData <- list()
        current.cohort <- cmPat/chSize
        if(cmPat < trialSize){
                cat("We are recommending the dose for your next cohort of patients...\n")
        }else{
                cat("The trial is ending and we are declaring efficacious doses, as well as recommending the lowest dose that is efficacious...\n")
        }
        cat(sprintf("The maximum sample size is : %d\nThe current enrolled number of patients are: %d\nThe current enrolled cohort is: %d\n", trialSize, cmPat, current.cohort))
        
        if(cmPat <= trialSize/2){
                if(cmPat <  trialSize/2){
                        cat("You are right now in stage 1, doing dose-escaltion based on safety only\n")
                }
                icycle <- with(toxicity.dat, ifelse(cycle == 1, 0, 1))
                PatData$toxicity <- data.frame(cbind(toxicity.dat[, c("DLT", "nTTP")], 1, 
                                                     toxicity.dat[, c("dose","cycle")], 
                                                     icycle))
                PatData$toxicity$subID <- toxicity.dat[, "subID"]
                names(PatData$toxicity) <- c("DLT", "nTTP", "intercept", "dose", "cycle", "icycle", "subID")
                ###################################################################################
                ################Model Estimation given PatData#####################################
                ################Stage 1 only considers the toxicity model##########################
                ###################################################################################
                nTTP <- PatData$toxicity[, 2]
                X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
                n.subj <- length(unique(PatData$toxicity[, "subID"]))
                W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                        which(a == unique(PatData$toxicity$subID))
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
                        ####draw posterior plots####
                        nTTP_post1 <-  sapply(doses, function(a){
                                apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a)
                                })
                        })
                        
                        nTTP1_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post1) / length(doses)), 
                                                Post_nTTP = c(nTTP_post1))
                        p_nTTP1 <- ggplot(data = nTTP1_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                                stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Cycle 1") + 
                                geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                                scale_x_discrete(breaks = doses, 
                                                 labels = doses)
                        
                        #----------------------------------------------------------------#
                        nTTP_post2 <-  sapply(doses, function(a){
                                apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3])
                                })
                        })
                        
                        nTTP2_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post2) / length(doses)), 
                                                Post_nTTP = c(nTTP_post2))
                        p_nTTP2 <- ggplot(data = nTTP2_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                                stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Late Cycles") + 
                                geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                                scale_x_discrete(breaks = doses, 
                                                 labels = doses)
                        #----------------------------------------------------------------#
                        
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
                        
                        ####toxicity profile 1######
                        toxpf1 <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a)
                                }))
                        })
                        
                        ####toxicity profile 2######
                        toxpf2 <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3] * 1)
                                }))
                        })
                        
                        tox.pf <- data.frame(rbind(toxpf1, toxpf2))
                        names(tox.pf) <- paste0("dose", doses)
                        row.names(tox.pf) <- c("toxicity_profile1", "toxicity_profile2")
                        
                        if(length(allow.doses) == 0){
                                cat("No doses are safe so the trial is terminated.\nThe toxicity profile for each doses are shown below:\n")
                                print(tox.pf)
                                return(results = list(nxtdose = 0, 
                                                      tox.pf = tox.pf, 
                                                      allow.doses = 0, 
                                                      p_nTTP1 = p_nTTP1,
                                                      p_nTTP2 = p_nTTP2))
                        }
                        cat("This is the end of stage 1, ready to enter stage 2.\nThe set of allowable (safe) doses based on safety only is as below:\n")
                        print(allow.doses)
                        cat("Posterior estimates (mean) of toxicity:\n")
                        print(tox.pf)
                        
                }else{
                        sim.betas <- as.matrix(rbind(post.samples$beta_other[1,,1], post.samples$beta_dose[,,1], 
                                                     post.samples$beta_other[2,,1]))
                        ####draw posterior plots####
                        nTTP_post1 <-  sapply(doses, function(a){
                                apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a)
                                })
                        })
                        
                        nTTP1_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post1) / length(doses)), 
                                                Post_nTTP = c(nTTP_post1))
                        p_nTTP1 <- ggplot(data = nTTP1_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                                stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Cycle 1") + 
                                geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                                scale_x_discrete(breaks = doses, 
                                                 labels = doses)
                        
                        #----------------------------------------------------------------#
                        nTTP_post2 <-  sapply(doses, function(a){
                                apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3])
                                })
                        })
                        
                        nTTP2_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post2) / length(doses)), 
                                                Post_nTTP = c(nTTP_post2))
                        p_nTTP2 <- ggplot(data = nTTP2_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                                stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Late Cycles") + 
                                geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                                scale_x_discrete(breaks = doses, 
                                                 labels = doses)
                        #----------------------------------------------------------------#
                        loss.doses <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        abs(o[1] + o[2] * a - tox.target)
                                }))
                        })
                        nxtdose <- doses[which.min(loss.doses)]
                        doseA <- with(toxicity.dat, dose[grepl(paste0("cohort", current.cohort, "subject"), subID)])[1]
                        
                        if(as.numeric(nxtdose) > (doseA + 1)){
                                nxtdose <- doseA + 1;
                        }
                        
                        ####toxicity profile 1######
                        toxpf1 <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a)
                                }))
                        })
                        
                        ####toxicity profile 2######
                        toxpf2 <- sapply(doses, function(a){
                                mean(apply(sim.betas, 2, function(o){
                                        as.numeric(o[1] + o[2] * a + o[3] * 1)
                                }))
                        })
                        
                        tox.pf <- data.frame(rbind(toxpf1, toxpf2))
                        names(tox.pf) <- paste0("dose", doses)
                        row.names(tox.pf) <- c("toxicity_profile1", "toxicity_profile2")
                        
                        if (cmPat == chSize){
                                DLT <- with(toxicity.dat, DLT[cycle == 1 & grepl(paste0("cohort", 1, "subject"), subID)])
                                if (sum(DLT) == 0){
                                        nxtdose <- 2
                                        dose_flag <- 0
                                }else{
                                        nxtdose <- 1
                                        dose_flag <- 1			
                                }
                        }
                        if ((cmPat >= chSize*2) & (dose_flag == 1)){
                                DLT <-  with(toxicity.dat, DLT[cycle == 1 & grepl(paste0("cohort", current.cohort, "subject"), subID)])
                                if (sum(DLT) == 0){
                                        nxtdose <- 2
                                        dose_flag <- 0
                                }else{
                                        nxtdose <- 1
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
                                        cat("No doses are safe so the trial is terminated.\nPosterior estimates (mean) of toxicity:\n")
                                        print(tox.pf)
                                        return(results = list(nxtdose = 0, 
                                                              tox.pf = tox.pf, 
                                                              allow.doses = 0,
                                                              p_nTTP1 = p_nTTP1, 
                                                              p_nTTP2 = p_nTTP2))
                                }
                        }
                        cat("The trial shall continue and posterior estimates (mean) of toxicity:\n")
                        print(tox.pf)
                        cat(sprintf("Next recommended dose: %d\n", nxtdose))
                        return(results = list(nxtdose = nxtdose, 
                                              tox.pf = tox.pf, 
                                              dose_flag = dose_flag,
                                              p_nTTP1 = p_nTTP1, 
                                              p_nTTP2 = p_nTTP2))
                }
                
        }
        
        if(cmPat >= trialSize/2 & cmPat < trialSize){
                if(cmPat == trialSize/2){
                        cat("You completed stage 1 and are entering stage 2..., randomizing the next cohort of patients towards higher predicted efficacy...\n")
                }else{
                        cat("You are right now in stage 2, randomizing the next cohort of patients towards higher predicted efficacy...\n")
                }
                icycle <- with(toxicity.dat, ifelse(cycle == 1, 0, 1))
                PatData$toxicity <- data.frame(cbind(toxicity.dat[, c("DLT", "nTTP")], 1, 
                                                     toxicity.dat[, c("dose","cycle")], 
                                                     icycle))
                PatData$toxicity$subID <- toxicity.dat[, "subID"]
                names(PatData$toxicity) <- c("DLT", "nTTP", "intercept", "dose", "cycle", "icycle", "subID")
                
                square.dose <- with(efficacy.dat, dose^2)
                PatData$efficacy <- data.frame(cbind(efficacy.dat[, "Efficacy"], 1, 
                                                     efficacy.dat[, "dose"], square.dose))
                PatData$efficacy$subID <- efficacy.dat[, "subID"]
                names(PatData$efficacy) <- c("Efficacy", "intercept", "dose", "square.dose", "subID")
                
                ############Model Estimation to randomize among allowable doses##########
                nTTP <- PatData$toxicity[, 2]
                X_y <- as.matrix(PatData$toxicity[, c("intercept", "dose", "icycle")])
                n.subj <- length(unique(PatData$toxicity[, "subID"]))
                W_y <- matrix(0, nrow(PatData$toxicity), n.subj)
                W_y[cbind(1:nrow(W_y), as.numeric(sapply(PatData$toxicity$subID, function(a){
                        which(a == unique(PatData$toxicity$subID))
                })))] <- 1
                
                EFF <- PatData$efficacy[, 1]
                X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "square.dose")])
                keepeff.ind <- sapply(PatData$efficacy$subID, function(a){
                        which(as.character(a) == unique(PatData$toxicity$subID))
                })
                
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
                
                ####draw posterior plots####
                nTTP_post1 <-  sapply(doses, function(a){
                        apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a)
                        })
                })
                
                nTTP1_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post1) / length(doses)), 
                                        Post_nTTP = c(nTTP_post1))
                p_nTTP1 <- ggplot(data = nTTP1_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Cycle 1") + 
                        geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                        scale_x_discrete(breaks = doses, 
                                         labels = doses)
                
                #----------------------------------------------------------------#
                nTTP_post2 <-  sapply(doses, function(a){
                        apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3])
                        })
                })
                
                nTTP2_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post2) / length(doses)), 
                                        Post_nTTP = c(nTTP_post2))
                p_nTTP2 <- ggplot(data = nTTP2_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Late Cycles") + 
                        geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                        scale_x_discrete(breaks = doses, 
                                         labels = doses)
                #----------------------------------------------------------------#
                sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                EFF_post <- sapply(doses, function(a){
                        apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        })
                })
                EFFp_dat <- data.frame(Dose = rep(doses, each = length(EFF_post) / length(doses)), 
                                       Post_EFF = c(EFF_post))
                p_EFF <- ggplot(data = EFFp_dat, aes(x = factor(Dose), y = Post_EFF, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("Posterior Efficacy") + 
                        xlab("Dose") + scale_x_discrete(breaks = doses, 
                                                        labels = doses)
                
                
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
                
                ####toxicity profile 1######
                toxpf1 <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a)
                        }))
                })
                
                ####toxicity profile 2######
                toxpf2 <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * 1)
                        }))
                })
                
                tox.pf <- data.frame(rbind(toxpf1, toxpf2))
                names(tox.pf) <- paste0("dose", doses)
                row.names(tox.pf) <- c("toxicity_profile1", "toxicity_profile2")
                ####efficacy profile####
                #sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                eff.pf <- sapply(doses, function(a){
                        mean(apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        }))
                })
                eff.pf <- data.frame(rbind(eff.pf))
                names(eff.pf) <- paste0("dose", doses)
                row.names(eff.pf) <- c("efficacy_profile")
                
                if(length(allow.doses) == 0){
                        if(cmPat == trialSize/2){
                                cat("No doses are safe so the trial is terminated.\n The updated posterior estimates (mean) of toxicity and efficacy:\n")
                        }else{
                                cat("No doses are safe so the trial is terminated.\n Posterior estimates (mean) of toxicity and efficacy:\n")
                        }
                        #cat("Posterior estimates (mean) of toxicity and efficacy:\n")
                        print(rbind(tox.pf, eff.pf))
                        return(results = list(nxtdose = 0, 
                                              tox.pf = tox.pf, 
                                              eff.pf = eff.pf,
                                              allow.doses = 0, 
                                              p_nTTP1 = p_nTTP1, 
                                              p_nTTP2 = p_nTTP2,
                                              p_EFF = p_EFF))
                }
                
                RAND.EFF <- sapply(allow.doses, function(a){
                        mean(apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        }))
                })
                RAND.EFF <- exp(RAND.EFF)/sum(exp(RAND.EFF))
                nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)
                
                ####a condition for untried higher dose level that is randomized#### 
                maxal_dose <- max(with(toxicity.dat, dose))
                if(nxtdose > maxal_dose + 1){
                        nxtdose <- maxal_dose + 1
                }
                if(cmPat == trialSize/2){
                        cat("The updated posterior estimates (mean) of toxicity and efficacy:\n")
                }else{
                        cat("Posterior estimates (mean) of toxicity and efficacy:\n")
                }
                #cat("Posterior estimates (mean) of toxicity and efficacy:\n")
                print(rbind(tox.pf, eff.pf))
                cat(sprintf("Next recommended dose: %d\n", nxtdose))
                return(results = list(nxtdose = nxtdose, 
                                      tox.pf = tox.pf, 
                                      eff.pf = eff.pf,
                                      allow.doses = allow.doses, 
                                      p_nTTP1 = p_nTTP1, 
                                      p_nTTP2 = p_nTTP2,
                                      p_EFF = p_EFF))
                
        }else if(cmPat == trialSize){
                cat("The trial is completed.\nYou completed stage 2 and are entering stage 3...\n")
                icycle <- with(toxicity.dat, ifelse(cycle == 1, 0, 1))
                PatData$toxicity <- data.frame(cbind(toxicity.dat[, c("DLT", "nTTP")], 1, 
                                                     toxicity.dat[, c("dose","cycle")], 
                                                     icycle))
                PatData$toxicity$subID <- toxicity.dat[, "subID"]
                names(PatData$toxicity) <- c("DLT", "nTTP", "intercept", "dose", "cycle", "icycle", "subID")
                
                square.dose <- with(efficacy.dat, dose^2)
                PatData$efficacy <- data.frame(cbind(efficacy.dat[, "Efficacy"], 1, 
                                                     efficacy.dat[, "dose"], square.dose))
                PatData$efficacy$subId <- efficacy.dat[, "subID"]
                names(PatData$efficacy) <- c("Efficacy", "intercept", "dose", "square.dose", "subID")
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
                X_e <- as.matrix(PatData$efficacy[, c("intercept", "dose", "square.dose")])
                keepeff.ind <- sapply(PatData$efficacy$subID, function(a){
                        which(as.character(a) == unique(PatData$toxicity$subID))
                })
                
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
                
                #sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                
                ####draw posterior plots####
                nTTP_post1 <-  sapply(doses, function(a){
                        apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a)
                        })
                })
                
                nTTP1_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post1) / length(doses)), 
                                        Post_nTTP = c(nTTP_post1))
                p_nTTP1 <- ggplot(data = nTTP1_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Cycle 1") + 
                        geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                        scale_x_discrete(breaks = doses, 
                                         labels = doses)
                
                #----------------------------------------------------------------#
                nTTP_post2 <-  sapply(doses, function(a){
                        apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3])
                        })
                })
                
                nTTP2_dat <- data.frame(Dose = rep(doses, each = length(nTTP_post2) / length(doses)), 
                                        Post_nTTP = c(nTTP_post2))
                p_nTTP2 <- ggplot(data = nTTP2_dat, aes(x = factor(Dose), y = Post_nTTP, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("nTTP Posterior Late Cycles") + 
                        geom_hline(yintercept = tox.target, color = "blue", linetype = 2) + xlab("Dose") + 
                        scale_x_discrete(breaks = doses, 
                                         labels = doses)
                #----------------------------------------------------------------#
                sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                EFF_post <- sapply(doses, function(a){
                        apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        })
                })
                EFFp_dat <- data.frame(Dose = rep(doses, each = length(EFF_post) / length(doses)), 
                                       Post_EFF = c(EFF_post))
                p_EFF <- ggplot(data = EFFp_dat, aes(x = factor(Dose), y = Post_EFF, group = factor(Dose))) + 
                        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab("Posterior Efficacy") + 
                        xlab("Dose") + scale_x_discrete(breaks = doses, 
                                                        labels = doses)
                
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
                                as.numeric(mean(sapply(2:MaxCycle, function(m){
                                        as.numeric(o[1] + o[2] * a + o[3] * m)
                                })) <= thrd2)
                        }))
                })
                
                allow.doses <- which(prob1.doses >= p1 & prob2.doses >= p2)
                
                ####toxicity profile 1######
                toxpf1 <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * 1)
                        }))
                })
                
                ####toxicity profile 2######
                toxpf2 <- sapply(doses, function(a){
                        mean(apply(sim.betas, 2, function(o){
                                as.numeric(mean(sapply(2:MaxCycle, function(m){
                                        as.numeric(o[1] + o[2] * a + o[3] * m)
                                })))
                        }))
                })
                
                tox.pf <- data.frame(rbind(toxpf1, toxpf2))
                names(tox.pf) <- paste0("dose", doses)
                row.names(tox.pf) <- c("toxicity_profile1", "toxicity_profile2")
                ####efficacy profile####
                #sim.alphas <- as.matrix(rbind(post.samples$alpha[,,1]))
                eff.pf <- sapply(doses, function(a){
                        mean(apply(sim.alphas, 2, function(o){
                                as.numeric(o[1] + o[2] * a + o[3] * a^2)
                        }))
                })
                eff.pf <- data.frame(rbind(eff.pf))
                names(eff.pf) <- paste0("dose", doses)
                row.names(eff.pf) <- c("efficacy_profile1")
                if(length(allow.doses) == 0){
                        cat("No doses are safe so the trial is terminated.\n")
                        cat("Posterior estimates (mean) of toxicity and efficacy:\n")
                        print(rbind(tox.pf, eff.pf))
                        return(results = list(opt.dose = 0, 
                                              tox.pf = tox.pf, 
                                              eff.pf = eff.pf,
                                              allow.doses = 0, 
                                              p_nTTP1 = p_nTTP1, 
                                              p_nTTP2 = p_nTTP2,
                                              p_EFF = p_EFF))
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
                cat("Posterior estimates (mean) of toxicity and efficacy:\n")
                print(rbind(tox.pf, eff.pf))
                cat("The efficacious doses are:\n")
                print(proxy.eff.doses)
                cat(sprintf("We recommend the lowest efficacious dose:%d\n", recom.doses))
                return(results = list(opt.dose = recom.doses, 
                                      tox.pf = tox.pf, 
                                      eff.pf = eff.pf,
                                      allow.doses = allow.doses, 
                                      p_nTTP1 = p_nTTP1, 
                                      p_nTTP2 = p_nTTP2,
                                      p_EFF = p_EFF))
                
        }
        
}

