
#convert the grade to DLT
#DLT: type 1 grade >=3. Or type 2 grade >=3. Or type 3 grade >=4.
grade2dlt <- function (w1_index,w2_index,w3_index){
dlt=rep(0,length(w1_index));
#because the first grade counts as grade 0, the criteria are 4,4,5 instead of 3,3,4
dlt[(w1_index >= 4) | (w2_index >= 4) | (w3_index >= 5)]=1;
return(dlt);
}

#calculate the average nTTP, and pr(DLT) of each dosage k
nTTPbar <- function (prob, v=2.5, sDose=1:6,
                     wm = matrix(c(0, 0.5, 0.75, 1, 1.5,
                                   0, 0.5, 0.75, 1, 1.5,
                                   0, 0, 0, 0.5, 1), byrow=TRUE, nrow=3, ncol=5)){
  nTTP=rep(NA, length(sDose));
  pDLT=rep(NA, length(sDose));
  for (k in sDose) {
    #probability of the grades
    #the outer product follow the order that for the first elements of the type 1, it goes through all elements of 2,3
    p=c(outer(prob$probaT3[k,],c(outer(prob$probaT2[k,],prob$probaT1[k,]))));
    #TTP of the grades
    w1<-wm[1,];w2<-wm[2,];w3<-wm[3,];
    w1_index=c(outer(rep(1,length(w3)),c(outer(rep(1,length(w2)),1:length(w1)))));
    this_w1=w1[w1_index];
    w2_index=c(outer(rep(1,length(w3)),c(outer(1:length(w2),rep(1,length(w1))))));
    this_w2=w2[w2_index];
    w3_index=c(outer(1:length(w3),c(outer(rep(1,length(w2)),rep(1,length(w1))))));
    this_w3=w3[w3_index];
    this_w=sqrt(this_w1^2+this_w2^2+this_w3^2);
    #dlt of the grades
    this_dlt=grade2dlt(w1_index,w2_index,w3_index);

    nTTP[k] = sum(p*this_w/v)
    pDLT[k] = sum(p*this_dlt)
  }
return(list(nTTP=nTTP, pDLT=pDLT))
}
