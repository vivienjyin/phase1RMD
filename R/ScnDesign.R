
#source('scenario_nTTB_nonNormalise_CDR4sc.R')
#source('ToxGrade2nTTP.R')
#source("transform_normal_grade-nTTB.R")

# toxicity probabilities for subsequent cycles;
# shift normal mean from 0 to a;
# given cut-offs on the normal pdf;
TP <- function(cut, shift) {
  n.dose = dim(cut)[1]
  max.grade = dim(cut)[2] # AE grade 0 to max.grade
  n.grade = max.grade + 1
  n.type = dim(cut)[3]
  
  proba <- array(data=NA, dim=c(n.dose, max.grade+1, n.type))
  for (d in 1:n.dose){ # sDose = 1:6
    for (t in 1:n.type) { # 3 types of AEs
      cumPr <- pnorm(cut[d,,t], mean=shift)
      celPr <- c(cumPr[1], diff(cumPr, lag=1, differences=1), 1-cumPr[max.grade])
      proba[d,,t] = celPr
    }
  }
  
  #prob <- list(probaT1=proba[,,1], probaT2=proba[,,2], probaT3=proba[,,3])
  return(proba)
}

## assume normal(0,1), find out the cut-offs;
cutoff <- function(prob){
  n.dose <- dim(prob$probaT1)[1]
  n.grade <- dim(prob$probaT1)[2]
  n.type = 3; # this may need changes depending on how many AE types we have;
  
  proba <- array(data=NA, dim=c(n.dose, n.grade, n.type)) # 6 dose levels, 5 grades, 3 AE types;
  proba[,,1]=prob$probaT1
  proba[,,2]=prob$probaT2
  proba[,,3]=prob$probaT3
  cut <- array(data=NA, dim=c(n.dose, n.grade-1, n.type))  # 6 dose levels, 4 cut-offs to generate 5 grades, 3 AE types;
  for(d in 1:n.dose){
    for(t in 1:n.type){
      pg = proba[d,,t]
      cut[d,,t] = c(qnorm(pg[1]), qnorm(pg[1]+pg[2]), qnorm(pg[1]+pg[2]+pg[3]), qnorm(1-pg[5])) 
      ## this may need changes depends on how many AE grades we have;
   }
  }
  return(cut)
}

#

# generate tox prob matrix for multiple cycles
GenScn <- function(sc, sc.title="Scenario F, TRD = 4", n.cycle = 6, trend = 0.1, sDose=1:6){

  prob = desc_sc(sc)
  cut = cutoff(prob)
  n.dose <- dim(prob$probaT1)[1]
  n.grade <- dim(prob$probaT1)[2]
  n.type = dim(cut)[3];
 
  proba.cycle <- array(data=NA, dim=c(n.dose, n.grade, n.type, n.cycle))
  proba.cycle[,,1,1] = prob$probaT1
  proba.cycle[,,2,1] = prob$probaT2
  proba.cycle[,,3,1] = prob$probaT3
  
  res = nTTPbar(prob=prob)
  #plot(x=sDose, y=as.vector(res$nTTP), xlab="Standardized dose", ylab="mean nTTP", type ="b", ylim=c(0,1), main = sc.title)
  #abline(h=0.28, col=1, lty=3)
  
  count.cycle = 1
  for (a in seq(trend, by=trend, length=(n.cycle-1))) { # a in seq(0,1, by=0.1, length=5);
    count.cycle = count.cycle + 1
    proba.late = TP(cut, a) # shift mean to from 0 to a in the normal distribution to generate subsequent cycle tox prob;
    proba.cycle[,,,count.cycle]=proba.late
    probc <- list(probaT1=proba.late[,,1], probaT2=proba.late[,,2], probaT3=proba.late[,,3])
    resc = nTTPbar(prob=probc)
    #print(res4c$nTTP)
    lines(x=sDose, y=as.vector(resc$nTTP), type ="b", col=2)
  }
  
  return(proba.cycle)
}

# generate tox prob matrix for multiple cycles
# Then simulate the patient nTTP based on the prob matrix
GenScn_nTTP_inside <- function(sc, sdose=1:6, n.cycle = 6, trend = 0.1){

  prob = desc_sc(sc)
  cut = cutoff(prob)
  n.dose <- dim(prob$probaT1)[1]
  n.grade <- dim(prob$probaT1)[2]
  n.type = dim(cut)[3];
  
  proba.cycle <- array(data=NA, dim=c(n.dose, n.grade, n.type, n.cycle))
  proba.cycle[,,1,1] = prob$probaT1
  proba.cycle[,,2,1] = prob$probaT2
  proba.cycle[,,3,1] = prob$probaT3
  
  #the nTTP of dose and cycle
  nTTP_res <- matrix(0,nrow=max(sdose),ncol=n.cycle)
  
  res = nTTPbar(prob=prob)
  #res = nTTP_dist(prob=prob, dose=dose)
  nTTP_res[,1]=res$nTTP;
  
  count.cycle = 1
  for (a in seq(trend, by=trend, length=(n.cycle-1))) { # a in seq(0,1, by=0.1, length=5);
    count.cycle = count.cycle + 1
    proba.late = TP(cut, a) # shift mean to from 0 to a in the normal distribution to generate subsequent cycle tox prob;
    proba.cycle[,,,count.cycle]=proba.late
    probc <- list(probaT1=proba.late[,,1], probaT2=proba.late[,,2], probaT3=proba.late[,,3])
    resc = nTTPbar(prob=probc)
	#resc = nTTP_dist(prob=probc, dose=dose)
    #print(res4c$nTTP)
    #lines(x=sDose, y=as.vector(resc$nTTP), type ="b", col=2)
	#nTTP_res=c(nTTP_res,resc);
	nTTP_res[,count.cycle]=resc$nTTP
  }

  return(nTTP_res)
}

GenScn_nTTP_probc <- function(sc, count.cycle = 1, trend = 0.1){

  prob = desc_sc(sc)
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

# generate tox prob matrix for multiple cycles
# Then simulate the patient nTTP based on the prob matrix
GenScn_nTTP <- function(sc, dose=1, count.cycle = 1, trend = 0.1){

  prob = desc_sc(sc)
  #print('prob');print(prob);
  cut = cutoff(prob)
  n.dose <- dim(prob$probaT1)[1]
  n.grade <- dim(prob$probaT1)[2]
  n.type = dim(cut)[3];
 
  #the nTTP of dose and cycle
  nTTP_res <- NULL;
  
  #res = nTTPbar(prob=prob)
  res = nTTP_dist(prob=prob, dose=dose)
  
  if (count.cycle == 1){
    return(res);
  }else{
    a <- trend * (count.cycle - 1);
    proba.late = TP(cut, a) # shift mean to from 0 to a in the normal distribution to generate subsequent cycle tox prob;
    probc <- list(probaT1=proba.late[,,1], probaT2=proba.late[,,2], probaT3=proba.late[,,3])
	resc = nTTP_dist(prob=probc, dose=dose)
  }
  #print('probc');print(probc);
  return(resc)
}

#convert the grade to DLT
#DLT: type 1 grade >=3. Or type 2 grade >=3. Or type 3 grade >=4.
grade2dlt <- function (w1_index,w2_index,w3_index){
dlt=rep(0,length(w1_index));
#because the first grade counts as grade 0, the criteria are 4,4,5 instead of 3,3,4
dlt[(w1_index >= 4) | (w2_index >= 4) | (w3_index >= 5)]=1;
return(dlt);
}

nTTP_dist <- function (prob, dose=1, v=2.5, 
                     wm = matrix(c(0, 0.5, 0.75, 1, 1.5,
                                   0, 0.5, 0.75, 1, 1.5,
                                   0, 0, 0, 0.5, 1), byrow=TRUE, nrow=3, ncol=5)){
								   
#just simulate one patient
N<-1;								   
#TTP of the grades
w1<-wm[1,];w2<-wm[2,];w3<-wm[3,];

#sample with for loop
this_w1<-NULL;this_w2<-NULL;this_w3<-NULL;
i<-1;
while (i <= N){
w1_N<-sample(1:5,size=1,prob=prob$probaT1[dose,],replace=TRUE)
w2_N<-sample(1:5,size=1,prob=prob$probaT2[dose,],replace=TRUE)
w3_N<-sample(1:5,size=1,prob=prob$probaT3[dose,],replace=TRUE)
#generate the TTP for each patient
this_w1<-c(this_w1,w1[w1_N]);
this_w2<-c(this_w2,w2[w2_N]);
this_w3<-c(this_w3,w3[w3_N]);
i=i+1;
}
dlt=0;dlt[(w1_N >= 4) | (w2_N >= 4) | (w3_N >= 5)]=1;
# ###########
this_ttp<-sqrt(this_w1^2+this_w2^2+this_w3^2);
this_nttp<-this_ttp/v;
return(list(nttp=this_nttp,dlt=dlt));
}

# ## An example to calculate the probability of scenariors 
# ##################################################################################################

# pdf("Scenarios_6cycles.pdf")

# par(mfrow=c(2,2));

# GenScn(6, sc.title="Scenario A, TRD = 2", trend = 0.1)
# GenScn(6, sc.title="Scenario A, TRD = 2", trend = -0.1)

# GenScn(5, sc.title="Scenario C, TRD = 3", trend = 0.1)
# GenScn(5, sc.title="Scenario C, TRD = 3", trend = -0.1)

# GenScn(4, sc.title="Scenario F, TRD = 4", trend = 0.1)
# GenScn(4, sc.title="Scenario F, TRD = 4", trend = -0.1)

# GenScn(7, sc.title="Scenario G, TRD = 5", trend = 0.1)
# GenScn(7, sc.title="Scenario G, TRD = 5", trend = -0.1)

# dev.off()

# ##################################################################################################

# pdf("Scenarios_1cycle.pdf")

# prob4=desc_sc(4);

# res4 = nTTPbar(prob=prob4)
# plot(x=sDose, y=as.vector(res4$nTTP), xlab="Standardized dose", ylab="mean nTTP", type ="b", ylim=c(0,1), 
     # main = "Scenarios, 1st Cycle", col=4)
# abline(h=0.28, col=1, lty=2)
# legend("topright", col=2:5, lty=1,
       # c("TRD=2, Scenario A", "TRD=3, Scenario C", "TRD=4, Scenario F", "TRD=5, Scenario G"))

# prob6=desc_sc(6);
# res6 = nTTPbar(prob6)
# lines(x=sDose, y=as.vector(res6$nTTP), type ="b", col=2)

# prob5=desc_sc(5);
# res5 = nTTPbar(prob5)
# lines(x=sDose, y=as.vector(res5$nTTP), type ="b", col=3)

# prob7=desc_sc(7);
# res7 = nTTPbar(prob7)
# lines(x=sDose, y=as.vector(res7$nTTP), type ="b", col=5)

# dev.off()

