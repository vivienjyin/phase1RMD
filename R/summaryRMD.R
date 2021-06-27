

summaryRMD <- function(fit, fit2, pathout, numTrials=1, thisTrial){

  dir.create(paste(pathout, "/", thisTrial, sep=""));
  cmPat <- fit$data$n.subj # number of subjects that has non-mising values
  
  # beta and sigma2;
  parms <- mcmc(cbind(fit$beta.dose, fit$beta.other, fit$s2.gamma, fit$s2.epsilon))
  if (ncol(fit$beta.other) == 2) {
    varnames(parms) <- c("beta.dose", "beta.intercept", "beta.cycle", "s2.gamma", "s2.epsilon")
  }else{
    varnames(parms) <- c("beta.dose", "beta.intercept", "s2.gamma", "s2.epsilon")
  }
  parms2 <- mcmc(cbind(fit2$beta.dose, fit2$beta.other, fit2$s2.gamma, fit2$s2.epsilon))
  if (ncol(fit2$beta.other) == 2) {
    varnames(parms2) <- c("beta.dose", "beta.intercept", "beta.cycle", "s2.gamma", "s2.epsilon")
  }else{
    varnames(parms2) <- c("beta.dose", "beta.intercept", "s2.gamma", "s2.epsilon")
  }
  combined=mcmc.list(mcmc(parms),mcmc(parms2));
  res = summary(parms)
  write.table(res$statistics, file=paste(pathout, "/", thisTrial, "/statsParms", cmPat, ".txt", sep=""))
  write.table(res$quantile, file=paste(pathout,  "/", thisTrial, "/quanParms", cmPat, ".txt", sep=""))
  write.table(effectiveSize(parms), file=paste(pathout,  "/", thisTrial, "/sizeParms", cmPat, ".txt", sep=""))

  # gamma
  gamma.mc <- mcmc(fit$gamma)
  colnames.gamma <- NULL
  for (i in 1:cmPat)
    colnames.gamma <- c(colnames.gamma, paste("gamma", i, sep=''))
  varnames(gamma.mc) <- colnames.gamma
  res = summary(gamma.mc)
  write.table(res$statistics, file=paste(pathout,  "/", thisTrial, "/statsGamma", cmPat, ".txt", sep=""))
  write.table(res$quantile, file=paste(pathout,  "/", thisTrial, "/quanGamma", cmPat, ".txt", sep=""))
  write.table(effectiveSize(gamma.mc), file=paste(pathout,  "/", thisTrial, "/sizeGamma", cmPat, ".txt", sep=""))

  #pdf(file=paste(pathout,  "/", thisTrial, "/traceParms", cmPat, ".pdf", sep=""))
  #plot(parms)
  #dev.off()  

  #pdf(file=paste(pathout,  "/", thisTrial, "/traceGamma", cmPat, ".pdf", sep=""))
  #plot(gamma.mc)
  #dev.off()
  
  #if (numTrials == 1){  
  if (thisTrial ==1 ){  #plot all
   # traceplots

   #pdf(file=paste(pathout,  "/", thisTrial, "/gelmanParms", cmPat, ".pdf", sep=""))
   #gelman.plot(combined)
   #dev.off();

  }
}
