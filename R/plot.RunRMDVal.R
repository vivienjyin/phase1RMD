#plot.RunRMDVal <- function(res,cycle = 1){
plot.RunRMDVal <- function(x, ...){
	thislist <- list(...); 
	if (length(thislist$cycle) == 0){cycle <- 1;
	}else{cycle <- thislist$cycle;}
	res <- x;
	#nTTP <- c(res$cycle1.nTTP) #cycle1.nTTP has dose X simulations
	cycle.index <- seq((cycle-1)*length(res$sdose)+1, (cycle)*length(res$sdose), 1);
	nTTP <- res$nTTP.p[cycle.index,];
	dose <- rep(res$sdose, times = dim(nTTP)[2])
	#boxplot(nTTP ~ dose, xlab='Dose',ylab='nTTP Posterior')
	#abline(h = 0.28)
	nTTP <- c(nTTP)
	dat <- data.frame(nTTP = nTTP, dose = dose);
	ggplot(data = dat, aes(x = factor(dose), y = nTTP, group = factor(dose))) + 
        stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() + ylab(paste("nTTP Posterior Cycle ",cycle,sep='')) + 
        geom_hline(yintercept = 0.28, color = "blue", linetype = 2) + xlab("Dose") + 
        scale_x_discrete(breaks = 1:6, labels = 1:6)	

}

