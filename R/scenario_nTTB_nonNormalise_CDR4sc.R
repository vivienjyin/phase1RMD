#rm(list=ls())
#setwd("D:/These/simulation_R/multiple_scenarios/dose-allocation algorithms/3algo/SimulationsLancer/algo-Rui")


###########################
# Definition of the weights

H      <- 3            # Type of toxicity h
G      <- 5            # Grade of toxicity
K      <- 6            # Dose j

listdo <- c("D1","D2","D3","D4","D5","D6")
listox <- c("Renal","Neuro","Hemato")
listgr <- c("G0","G1","G2","G3","G4")
listt1 <- c("R0","R1","R2","R3","R4")
listt2 <- c("N0","N1","N2","N3","N4")
listt3 <- c("H0","H1","H2","H3","H4")

dimarr <- list(listt1,listt2,listt3)

# Array of the "name" of the 125 combinations of toxicities 
#----------------------------------------------------------
mtox   <- array(NA,c(5,5,5))

for (t1 in 1:5)
  {
  for (t2 in 1:5)
    {
    for (t3 in 1:5)
      {
      mtox[t1,t2,t3] <-  paste(listt1[t1],"_",listt2[t2],"_",listt3[t3],sep="")
      }
    }
  }

#--------------------------------------------------
# Matrix of weights for each grade/type of toxicity
#--------------------------------------------------

wm   <- matrix(NA,H,G)      
rownames(wm) <- listox
colnames(wm) <- listgr

#------------------------------------------------
# DATA ENTRY: Definition of the matrix of weights
#------------------------------------------------
w1     <- c(0, 0.5, 0.75, 1  , 1.5)  #Renal
w2     <- c(0, 0.5, 0.75, 1  , 1.5)  #Neuro 
w3     <- c(0, 0  , 0   , 0.5, 1  )  #Hemato

#toxmax <- 2.5          # Value used to normalize the nttb 
toxmax  <- 1            # Value used to not normalize the nttb 
#------------------------------------------------
# End DATA ENTRY matrix of weights
#------------------------------------------------
wm <- rbind(w1,w2,w3)

# Array of the normalized scores associated with the 125 combinations 
#--------------------------------------------------------------------

###############################################################
# Vector of proba of grades 0, 1, 2, 3 & 4 for each type of tox

# For each dose and each tox type, 
# the proba of grade 4 is derived (1-sum(proba of grades(0:3))
proba_gr4 <- function(j,h)
    {
    case  <- paste("proba",j,h,sep="")
    proba <- get(case)
    return(c(proba,1-sum(proba))) 
    }

# Proba of DLT = Proba of [Renal gr(3,4) or Neuro gr(3,4) or Hem gr(4)]
pDLT     <-function(proba)
    {    
    return(sum(proba[c(4,5),      , ]) 
         + sum(proba[      ,c(4,5), ]) 
         + sum(proba[      ,      ,5])
         - sum(proba[c(4,5),,])*sum(proba[,c(4,5),])
         - sum(proba[c(4,5),,])*sum(proba[,,5])
         - sum(proba[,c(4,5),])*sum(proba[,,5])
         + sum(proba[c(4,5),,])*sum(proba[,c(4,5),])*sum(proba[,,5]))
    }
#Otherwise    
pDLTbis  <-function(proba)
    {    
    return(sum(proba[c(4,5),     , ]) 
         + sum(proba[(1:3),c(4,5), ]) 
         + sum(proba[(1:3),(1:3) ,5]))
             
    }


#----------------------------------------------------
# Description of the SCENARIOS
#----------------------------------------------------



desc_sc <- function(sc)
{

nTTB  <- function(w1,w2,w3)
{
nTTB <- array(NA,c(5,5,5))
dimnames(nTTB) <- dimarr
for (t1 in 1:5)
  {
  for (t2 in 1:5)
    {
    for (t3 in 1:5)
      {
      nTTB[t1,t2,t3] <-  ((w1[t1]^2+w2[t2]^2+w3[t3]^2)^0.5)/toxmax
      }
    }
  }
return(nTTB)
}
nTTB <- nTTB(w1,w2,w3)

# Mean(nTTB)
mnTTB    <-function(proba)
    {    
    return(sum(proba*nTTB))
    }

#----------------------------------------------------
# DATA ENTRY: Definition of the proba in the scenario
#----------------------------------------------------
# Name of the vector of proba(j,h)
# 1st digit: dose j
# 2nd digit: tox h

#revision031917
# <- to <-
if (sc == 1) 
{
#Scenario 1 - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05) #pg_m1[60,(2:6)] 
    proba21 <- c(0.6700314463,0.2550348541,0.0675900686,0.0070527739,0.0002908571) #pg_m1[79,(2:6)]  
    proba31 <- c(0.415440537,0.363105270,0.187570075,0.030628442,0.003255677) #pg_m1[109,(2:6)]    
    proba41 <- c(0.2291283,0.410444821,0.284169688,0.060550194,0.01570703) #pg_m1[118,(2:6)]
    proba51 <- c(0.1974597,0.428924923,0.295730337,0.06197536,0.015909665) #pg_m1[126,(2:6)]
    proba61 <- c(0.1826522,0.40126360,0.28506666,0.05724707,0.07377044) #pg_m1[135,(2:6)]
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
 #Scenario 1 - For the 6 doses - Tox 2
    proba12 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)  #pg_m1[8,(2:6)] 
    proba22 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)  #pg_m1[60,(2:6)]
    proba32 <- c(0.5900114463,0.3050348541,0.0975900686,0.0070527739,0.0003108571)  #pg_m1[79,(2:6)]
    proba42 <- c(0.3234405,0.363105270,0.2780202,0.030628442,0.0048055676) #pg_m1[109,(2:6)]     
    proba52 <- c(0.2798375,0.382924923,0.291730337,0.040597536,0.004909665) #pg_m1[118,(2:6)]
    proba62 <- c(0.1486522,0.37726360,0.42506666,0.04324707,0.00577044)  #pg_m1[126,(2:6)]
   probaT2  <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 1 - For the 6 doses - Tox 3
    proba13 <- c(9.221962e-01,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-06)   #pg_m1[30,(2:6)]
    proba23 <- c(0.351972708,0.380398399,0.215012755,0.048219650,0.004396488)  #pg_m1[120,(2:6)]
    proba33 <- c(0.03992951,0.15872049,0.36291812,0.30459721,0.13383468)  #pg_m1[188,(2:6)]
    proba43 <- c(0.02562873,0.1033714,0.30957618,0.32531566,0.23610805) #pg_m1[210,(2:6)] 
    proba53 <- c(0.02494002,0.032113417,0.223903876,0.381799964,0.337242727)  #pg_m1[230,(2:6)]
    proba63 <- c(0.012441242,0.012410522,0.150691295,0.364435254,0.460021686)  #pg_m1[250,(2:6)]
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)

}    
if (sc == 2) 
{
#pg_m1 modified
#Scenario 2 - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.297832e-02,1.334457e-03,6.672584e-05)
    proba21 <- c(0.6700314463,0.2550348541,0.0715900686,0.0030527739,0.0002908571)
    proba31 <- c(0.587064423,0.301703140,0.101023054,0.009568431,0.000640953)
    proba41 <- c(0.436440537,0.363105270,0.177570075,0.020628442,0.002255677)    
    proba51 <- c(0.3385375,0.382924923,0.251730337,0.022597536,0.004209665)
    proba61 <- c(0.12825223,0.38726360,0.38106666,0.05324707,0.05017044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 2 - For the 6 doses - Tox 2
    proba12 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)
    proba22 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)
    proba32 <- c(0.7323711065,0.2150127550,0.0482196501,0.0042491868,0.0001473015)
    proba42 <- c(0.6700314463,0.2550348541,0.0675900686,0.0070527739,0.0002908571)
    proba52 <- c(0.446928264,0.388444821,0.096169688,0.064550194,0.003907033)
    proba62 <- c(0.228537539,0.402924923,0.271730337,0.070597536,0.026209665)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 2 - For the 6 doses - Tox 3
    proba13 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)#pg_m1[8,(2:6)]   
    proba23 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)#pg_m1[60,(2:6)] 
    proba33 <- c(0.6700314463,0.2550348541,0.0675900686,0.0070527739,0.0002908571)#pg_m1[79,(2:6)]  
    proba43 <- c(0.436440537,0.363105270,0.167570075,0.030628442,0.002255677)#pg_m1[109,(2:6)]  
    proba53 <- c(0.366928264,0.378444821,0.206169688,0.044550194,0.003907033)#pg_m1[118,(2:6)] 
    proba63 <- c(0.317537539,0.382924923,0.241730337,0.050597536,0.007209665)#pg_m1[126,(2:6)] 
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}

if (sc == 3) 
{
#pg_m1 modified
#Scenario 3 - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)
    proba21 <- c(0.7600314463,0.1850348541,0.0475900686,0.0070527739,0.0002908571)
    proba31 <- c(0.7084405, 0.1991053,0.0505701,0.040628442,0.001255677)
    proba41 <- c(0.6464405,0.233105270,0.067570075,0.047628442,0.005255677)    
    proba51 <- c(0.397537539,0.382924923,0.151730337,0.060597536,0.007209665)
    proba61 <- c(0.2602522,0.37726360,0.28106666,0.07334707,0.00807044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 2
    proba12 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)
    proba22 <- c(0.816087,1.717286e-01,0.0097832,2.334457e-03,6.672584e-05)
    proba32 <- c(0.7787061,0.183105270,0.0099701,0.0069628442,0.001255677)
    proba42 <- c(0.7590314463,0.2050348541,0.0275900686,0.0070527739,0.0012908571)
    proba52 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba62 <- c(0.2602522,0.37726360,0.28106666,0.07324707,0.00817044)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 3
    proba13 <- c(9.221962e-01,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-06)
    proba23 <- c(0.656972708,0.280398399,0.010012755,0.048219650,0.004396488)
    proba33 <- c(0.58092951,0.20872049,0.02291812,0.10359721,0.08383468)
    proba43 <- c(0.01462873,0.19437138,0.23957618,0.29531566,0.25610805)
    proba53 <- c(0.004940016,0.052113417,0.223903876,0.381799964,0.337242727)
    proba63 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    


if (sc == 33) 
{
#pg_m1 modified
#Scenario 3bis - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)
    proba21 <- c(0.7600314463,0.1850348541,0.0475900686,0.0070527739,0.0002908571)
    proba31 <- c(0.6787845, 0.2251053,0.0655701,0.02628442,0.004255677)
    proba41 <- c(0.6664405,0.233105270,0.067570075,0.027628442,0.005255677)    
    #proba51 <- c(0.267537539,0.452924923,0.201730337,0.070597536,0.007209665)
    proba51 <- c(0.2802522,0.37726360,0.24106666,0.09334707,0.00807044)
    proba61 <- c(0.1592522,0.38726360,0.3100667,0.13334707,0.01007044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 2
    proba12 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)
    proba22 <- c(0.816087,1.717286e-01,0.0097832,2.334457e-03,6.672584e-05)
    proba32 <- c(0.7677061,0.197105270,0.0269701,0.0069628442,0.001255677)
    proba42 <- c(0.7590314463,0.2050348541,0.0275900686,0.0070527739,0.0012908571)
    #proba52 <- c(0.277037539,0.367924923,0.276730337,0.070597536,0.007709665)
    proba52 <- c( 0.2292522,0.3502236,0.318106666,0.09324707,0.00917044)
    proba62 <- c(0.1842522,0.36726360,0.3296666,0.10724707,0.01157044)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 3
    proba13 <- c(9.221962e-01,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-06)
    proba23 <- c(0.656972708,0.280398399,0.010012755,0.048219650,0.004396488)
    proba33 <- c(0.25092951,0.26872049,0.18291812,0.14359721,0.15383468)
    proba43 <- c(0.01462873,0.19437138,0.23957618,0.29531566,0.25610805)
    proba53 <- c(0.004940016,0.032113417,0.223903876,0.351799964,0.387242727)
    proba63 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    


if (sc == 333) 
{
#pg_m1 modified
#Scenario 3bisbis - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)
    proba21 <- c(0.7600314463,0.1850348541,0.0475900686,0.0070527739,0.0002908571)
    proba31 <- c(0.6787845, 0.2251053,0.0655701,0.02628442,0.004255677)
    proba41 <- c(0.4964405,0.303105270,0.147570075,0.047628442,0.005255677)    
    #proba51 <- c(0.267537539,0.452924923,0.201730337,0.070597536,0.007209665)
    proba51 <- c(0.2802522,0.37726360,0.24106666,0.09334707,0.00807044)
    proba61 <- c(0.1592522,0.38726360,0.3100667,0.13334707,0.01007044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 2
    proba12 <- c(9.685572e-01,2.932456e-02,2.061512e-03,5.610658e-05,5.869288e-07)
    proba22 <- c(0.816087,1.717286e-01,0.0097832,2.334457e-03,6.672584e-05)
    proba32 <- c(0.7677061,0.197105270,0.0269701,0.0069628442,0.001255677)
    proba42 <- c(0.5030314463,0.2750348541,0.1875900686,0.0300527739,0.0042908571)
    #proba52 <- c(0.277037539,0.367924923,0.276730337,0.070597536,0.007709665)
    proba52 <- c( 0.2292522,0.3502236,0.318106666,0.09324707,0.00917044)
    proba62 <- c(0.1842522,0.36726360,0.3296666,0.10724707,0.01157044)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 3 - For the 6 doses - Tox 3
    proba13 <- c(9.221962e-01,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-06)
    proba23 <- c(0.656972708,0.280398399,0.010012755,0.048219650,0.004396488)
    proba33 <- c(0.25092951,0.26872049,0.18291812,0.14359721,0.15383468)
    proba43 <- c(0.01462873,0.19437138,0.23957618,0.29531566,0.25610805)
    proba53 <- c(0.004940016,0.032113417,0.223903876,0.351799964,0.387242727)
    proba63 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    

if (sc == 4) 
{
#pg_m1 modified
#Scenario 4 - For the 6 doses - Tox 1
    proba11 <- c(0.7912914,1.717286e-01,3.197832e-02,4.334457e-03,6.672584e-04)
    proba21 <- c(0.7578314463,0.1800348541,0.0430900686,0.0100527739,0.0089908571)
    proba31 <- c(0.6854405, 0.1901053,0.0675701,0.043628442,0.013255677)
    proba41 <- c(0.6624405,0.200105270,0.077570075,0.045628442,0.014255677)    
    proba51 <- c(0.604537539,0.222924923,0.081730337,0.070597536,0.020209665)
    proba61 <- c(0.3902522,0.30726360,0.20106666,0.07334707,0.02807044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 4 - For the 6 doses - Tox 2
    proba12 <- c(0.9676741,2.932456e-02,2.061512e-03,8.610658e-04,7.869288e-05)
    proba22 <- c(0.813057,1.717286e-01,0.0057832,9.334457e-03,9.672584e-05)
    proba32 <- c(0.7617061,0.183105270,0.0409701,0.0099628442,0.004255677)
    proba42 <- c(0.6710314463,0.2050348541,0.1084900686,0.0100527739,0.0053908571)
    proba52 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba62 <- c(0.2602522,0.37726360,0.28106666,0.07324707,0.00817044)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 4 - For the 6 doses - Tox 3
    proba13 <- c(0.9172662,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-03)
    proba23 <- c(0.651972708,0.280398399,0.010012755,0.020219650,0.037396488)
    proba33 <- c(0.53592951,0.20872049,0.03091812,0.09059721,0.13383468)
    proba43 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba53 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba63 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    

#translate sc=4 to the left,TRD=3 instead of 4
if (sc == 5) 
{
#pg_m1 modified
#Scenario 5 - For the 6 doses - Tox 1
    proba11 <- c(0.7578314463,0.1800348541,0.0430900686,0.0100527739,0.0089908571)
    proba21 <- c(0.6854405, 0.1901053,0.0675701,0.043628442,0.013255677)
    proba31 <- c(0.6624405,0.200105270,0.077570075,0.045628442,0.014255677)    
    proba41 <- c(0.604537539,0.222924923,0.081730337,0.070597536,0.020209665)
    proba51 <- c(0.3902522,0.30726360,0.20106666,0.07334707,0.02807044)
    proba61 <- c(0.2892522,0.35726360,0.22106666,0.10334707,0.02907044)    
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 2
    proba12 <- c(0.813057,1.717286e-01,0.0057832,9.334457e-03,9.672584e-05)
    proba22 <- c(0.7617061,0.183105270,0.0409701,0.0099628442,0.004255677)
    proba32 <- c(0.6710314463,0.2050348541,0.1084900686,0.0100527739,0.0053908571)
    proba42 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba52 <- c(0.2602522,0.37726360,0.28106666,0.07324707,0.00817044)
    proba62 <- c(0.2082522,0.38726360,0.30106666,0.09324707,0.01017044)    
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 3
    proba13 <- c(0.651972708,0.280398399,0.010012755,0.020219650,0.037396488)
    proba23 <- c(0.53592951,0.20872049,0.03091812,0.09059721,0.13383468)
    proba33 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba43 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba53 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    proba63 <- c(0.004341242,0.020410522, 0.2096913,0.334435254, 0.4311217)    
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    



if (sc == 55) 
{
#pg_m1 modified
#Scenario 5 - For the 6 doses - Tox 1
    proba11 <- c(0.7908314463,0.1800348541,0.0130900686,0.0070527739,0.0089908571)
    proba21 <- c(0.5834405, 0.2201053,0.0775701,0.103628442,0.015255677)
    proba31 <- c(0.5604405,0.240105270,0.087570075,0.095628442,0.016255677)    
    proba41 <- c(0.391537539,0.282924923,0.161730337,0.143597536,0.020209665)
    proba51 <- c(0.2372522,0.32726360,0.20106666,0.19634707,0.03807044)
    proba61 <- c(0.1492522,0.35726360,0.22106666,0.21334707,0.05907044)    
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 2
    proba12 <- c(0.818077,1.717286e-01,0.0057832,4.334457e-03,7.672584e-05)  
    proba22 <- c(0.7017061,0.183105270,0.1009701,0.0099628442,0.004255677)
    proba32 <- c(0.6110314463,0.2050348541,0.1484900686,0.0300527739,0.0053908571)
    proba42 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba52 <- c(0.2442522,0.37726360,0.28106666,0.08924707,0.00817044)
    proba62 <- c(0.1002522,0.40726360,0.32406666,0.14924707,0.01917044)    
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 3
    proba13 <- c(0.591972708,0.240398399,0.090012755,0.040219650,0.037396488)
    proba23 <- c(0.37592951,0.20872049,0.13091812,0.15059721,0.13383468)
    proba33 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba43 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba53 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    proba63 <- c(0.004341242,0.020410522, 0.2096913,0.334435254, 0.4311217)    
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    

if (sc == 555) 
{
#pg_m1 modified
#Scenario 5 - For the 6 doses - Tox 1
    proba11 <- c(0.7908314463,0.1800348541,0.0130900686,0.0070527739,0.0089908571)
    proba21 <- c(0.6590405, 0.2201053,0.0575701,0.053028442,0.010255677)
    proba31 <- c(0.5624405,0.240105270,0.087570075,0.095628442,0.014255677)    
    proba41 <- c(0.391537539,0.282924923,0.161730337,0.143597536,0.020209665)
    proba51 <- c(0.2372522,0.32726360,0.20106666,0.19634707,0.03807044)
    proba61 <- c(0.1492522,0.35726360,0.22106666,0.21334707,0.05907044)    
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 2
    proba12 <- c(0.818077,1.717286e-01,0.0057832,4.334457e-03,7.672584e-05)  
    proba22 <- c(0.7217061,0.173105270,0.0909701,0.0099628442,0.004255677)
    proba32 <- c(0.6110314463,0.2050348541,0.1484900686,0.0300527739,0.0053908571)
    proba42 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba52 <- c(0.2442522,0.37726360,0.28106666,0.08924707,0.00817044)
    proba62 <- c(0.1002522,0.40726360,0.32406666,0.14924707,0.01917044)    
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 5 - For the 6 doses - Tox 3
    proba13 <- c(0.591972708,0.240398399,0.090012755,0.040219650,0.037396488)
    proba23 <- c(0.39592951,0.20872049,0.13091812,0.15059721,0.11383468)
    proba33 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba43 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba53 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    proba63 <- c(0.004341242,0.020410522, 0.2096913,0.334435254, 0.4311217)    
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}

#translate sc=4 to the left,TRD=2 instead of 4
if (sc == 6) 
{
#pg_m1 modified
#Scenario 6 - For the 6 doses - Tox 1
    proba11 <- c(0.6854405, 0.1901053,0.0675701,0.043628442,0.013255677)
    proba21 <- c(0.6624405,0.200105270,0.077570075,0.045628442,0.014255677)    
    proba31 <- c(0.604537539,0.222924923,0.081730337,0.070597536,0.020209665)
    proba41 <- c(0.3902522,0.30726360,0.20106666,0.07334707,0.02807044)
    proba51 <- c(0.2892522,0.35726360,0.22106666,0.10334707,0.02907044)  
    proba61 <- c(0.2482522,0.36726360,0.24106666,0.11334707,0.03007044)          
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 6 - For the 6 doses - Tox 2
    proba12 <- c(0.7617061,0.183105270,0.0409701,0.0099628442,0.004255677)
    proba22 <- c(0.6710314463,0.2050348541,0.1084900686,0.0100527739,0.0053908571)
    proba32 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    proba42 <- c(0.2602522,0.37726360,0.28106666,0.07324707,0.00817044)
    proba52 <- c(0.2082522,0.38726360,0.30106666,0.09324707,0.01017044)  
    proba62 <- c(0.1712522,0.40726360,0.31106666,0.09824707,0.01217044)          
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 6 - For the 6 doses - Tox 3
    proba13 <- c(0.53592951,0.20872049,0.03091812,0.09059721,0.13383468)
    proba23 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba33 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba43 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    proba53 <- c(0.004341242,0.020410522, 0.2096913,0.334435254, 0.4311217)    
    proba63 <- c(0.004141242,0.021410522, 0.2294913,0.324435254,0.4205217)        
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}    

#translate sc=4 to the right,TRD=5 instead of 4
if (sc == 7) 
{
#pg_m1 modified
#Scenario 7 - For the 6 doses - Tox 1
    proba11 <- c(0.8233914,1.517286e-01,2.197832e-02,2.334457e-03,5.672584e-04)    
    proba21 <- c(0.7912914,1.717286e-01,3.197832e-02,4.334457e-03,6.672584e-04)
    proba31 <- c(0.7578314463,0.1800348541,0.0430900686,0.0100527739,0.0089908571)
    proba41 <- c(0.6854405, 0.1901053,0.0675701,0.043628442,0.013255677)
    proba51 <- c(0.6624405,0.200105270,0.077570075,0.045628442,0.014255677)    
    proba61 <- c(0.604537539,0.222924923,0.081730337,0.070597536,0.020209665)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 7 - For the 6 doses - Tox 2
    proba12 <- c(0.9704214,2.672456e-02,2.009512e-03,7.688658e-04,7.569288e-05)   
    proba22 <- c(0.9676741,2.932456e-02,2.061512e-03,8.610658e-04,7.869288e-05)
    proba32 <- c(0.813057,1.717286e-01,0.0057832,9.334457e-03,9.672584e-05)
    proba42 <- c(0.7617061,0.183105270,0.0409701,0.0099628442,0.004255677)
    proba52 <- c(0.6710314463,0.2050348541,0.1084900686,0.0100527739,0.0053908571)
    proba62 <- c(0.397037539,0.257924923,0.276730337,0.060597536,0.007709665)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 7 - For the 6 doses - Tox 3
    proba13 <- c(0.930074,6.004359e-02,5.447148e-03,3.001706e-04,4.135045e-03) 
    proba23 <- c(0.9172662,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-03)
    proba33 <- c(0.651972708,0.280398399,0.010012755,0.020219650,0.037396488)
    proba43 <- c(0.53592951,0.20872049,0.03091812,0.09059721,0.13383468)
    proba53 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba63 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
    }

if (sc == 8) 
{
#pg_m1 modified
#Scenario 8 - For the 6 doses - Tox 1
    proba11 <- c(7.938919e-01,1.717286e-01,3.197832e-02,2.334457e-03,6.672584e-05)
    proba21 <- c(0.7578314463,0.1850348541,0.0430900686,0.0100527739,0.0039908571)
    proba31 <- c(0.6854405, 0.1991053,0.0675701,0.043628442,0.004255677)
    proba41 <- c(0.6454405,0.203105270,0.097570075,0.047628442,0.006255677)    
    proba51 <- c(0.557537539,0.212924923,0.151730337,0.070597536,0.007209665)
    proba61 <- c(0.4902522,0.22726360,0.18106666,0.09334707,0.00807044)
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 8 - For the 6 doses - Tox 2
    proba12 <- c(9.685270e-01,2.932456e-02,2.061512e-03,8.610658e-05,7.869288e-07)
    proba22 <- c(0.813057,1.717286e-01,0.0057832,9.334457e-03,9.672584e-05)
    proba32 <- c(0.7617061,0.183105270,0.0409701,0.0099628442,0.004255677)
    proba42 <- c(0.6930314463,0.2050348541,0.0864900686,0.0100527739,0.0053908571)
    proba52 <- c(0.607037539,0.217924923,0.106730337,0.060597536,0.007709665)
    proba62 <- c(0.5502522,0.24726360,0.12106666,0.07324707,0.00817044)
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 8 - For the 6 doses - Tox 3
    proba13 <- c(9.221962e-01,7.004359e-02,7.447148e-03,3.081706e-04,4.935045e-06)
    proba23 <- c(0.651972708,0.280398399,0.010012755,0.020219650,0.037396488)
    proba33 <- c(0.53592951,0.20872049,0.03091812,0.09059721,0.13383468)
    proba43 <- c(0.01462873,0.13437138,0.23957618,0.33531566,0.27610805)
    proba53 <- c(0.004940016,0.052113417,0.223903876,0.371799964,0.347242727)
    proba63 <- c(0.004441242,0.022410522, 0.2196913,0.344435254,0.409021686)
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}

if (sc == 9) 
{
#pg_m1 modified
#Scenario 9 - For the 6 doses - Tox 1
    proba11 <- c(0.4602522,0.25726360,0.18106666,0.07324707,0.02817044)
    proba21 <- c(0.3080000,0.35700000,0.20110000,0.07520000,0.05870000)
    proba31 <- c(0.1605338,0.37636000,0.296106666,0.0770000,0.09)
    proba41 <- c(0.1205338,0.34636000,0.35610620,0.0870000,0.09)
    proba51 <- c(0.1005338,0.33636000,0.336106666,0.067,0.1599995)    
    proba61 <- c(0.0805338,0.32636000,0.316106666,0.060,0.2069995)        
    probaT1 <-  round( matrix(c(proba11,proba21,proba31,proba41,proba51,proba61),6,5,byrow=T),digits=3)
#Scenario 9 - For the 6 doses - Tox 2
    proba12 <- c(0.4302522,0.30726360,0.14106666,0.11324707,0.00817044)
    proba22 <- c(0.2325803,0.3926360,0.20106666,0.164,0.009717044)
    proba32 <- c(0.099306,0.3876360,0.296106666,0.18524707,0.0317044)
    proba42 <- c(0.0705338,0.29636000,0.35610620,0.1870000,0.09)
    proba52 <- c(0.0505338,0.27636000,0.316106666,0.167,0.1899995)    
    proba62 <- c(0.0495338,0.27136000,0.310106666,0.164, 0.2049995)    
    probaT2 <-  round( matrix(c(proba12,proba22,proba32,proba42,proba52,proba62),6,5,byrow=T),digits=3)
#Scenario 9 - For the 6 doses - Tox 3
    proba13 <- c(0.202873,0.09437138,0.131332000000,0.33531566,0.23610805)
    proba23 <- c(0.130000,0.10000000,0.140000000000,0.39000000,0.24)
    proba33 <- c(0.035220,0.14000000,0.150000000000,0.4000000,0.27478)
    proba43 <- c(0.003000,0.15000000,0.151000000000,0.4040000,0.292)
    proba53 <- c(0.000000,0.15100000,0.152000000000,0.405000,0.292)
    proba63 <- c(0.000000,0.14600000,0.153000000000,0.408000,0.293)    
    probaT3 <-  round( matrix(c(proba13,proba23,proba33,proba43,proba53,proba63),6,5,byrow=T),digits=3)
}
#Proba of each of the 125 combinations of tox for each dose
#----------------------------------------------------------
proba1  <- array(NA,c(5,5,5) )
proba2  <- array(NA,c(5,5,5) )
proba3  <- array(NA,c(5,5,5) )
proba4  <- array(NA,c(5,5,5) )
proba5  <- array(NA,c(5,5,5) )
proba6  <- array(NA,c(5,5,5) )

# Dose 1
#-------
#proba11 <- proba_gr4(1,1)[-6]
#proba12 <- proba_gr4(1,2)[-6]
#proba13 <- proba_gr4(1,3)[-6]
proba1  <- outer(outer(proba11,proba12),proba13)
dimnames(proba1) <- dimarr
mnTTB1  <- mnTTB(proba1)
pDLT1   <- pDLT(proba1)
#pDLT1bis<- pDLTbis(proba1)

# Dose 2
#-------
#proba21 <- proba_gr4(2,1)[-6]
#proba22 <- proba_gr4(2,2)[-6]
#proba23 <- proba_gr4(2,3)[-6]
proba2  <- outer(outer(proba21,proba22),proba23)
dimnames(proba2) <- dimarr
mnTTB2  <- mnTTB(proba2)
pDLT2   <- pDLT(proba2)
#pDLT2bis<- pDLTbis(proba2)

# Dose 3
#-------
#proba31 <- proba_gr4(3,1)[-6]
#proba32 <- proba_gr4(3,2)[-6]
#proba33 <- proba_gr4(3,3)[-6]
proba3  <- outer(outer(proba31,proba32),proba33)
dimnames(proba3) <- dimarr
mnTTB3  <- mnTTB(proba3)
pDLT3   <- pDLT(proba3)
#pDLT3bis<-  pDLTbis(proba3)

# Dose 4
#-------
#proba41 <- proba_gr4(4,1)[-6]
#proba42 <- proba_gr4(4,2)[-6]
#proba43 <- proba_gr4(4,3)[-6]
proba4  <- outer(outer(proba41,proba42),proba43)
dimnames(proba4) <- dimarr
mnTTB4  <- mnTTB(proba4)
pDLT4   <- pDLT(proba4)
#pDLT4bis<-  pDLTbis(proba4)

# Dose 5
#-------
#proba51 <- proba_gr4(5,1)[-6]
#proba52 <- proba_gr4(5,2)[-6]
#proba53 <- proba_gr4(5,3)[-6]
proba5  <- outer(outer(proba51,proba52),proba53)
dimnames(proba5) <- dimarr
mnTTB5  <- mnTTB(proba5)
pDLT5   <- pDLT(proba5)
#pDLT5bis<-  pDLTbis(proba5)

# Dose 6
#-------
#proba61 <- proba_gr4(6,1)[-6]
#proba62 <- proba_gr4(6,2)[-6]
#proba63 <- proba_gr4(6,3)[-6]
proba6  <- outer(outer(proba61,proba62),proba63)
dimnames(proba6) <- dimarr
mnTTB6  <- mnTTB(proba6)
pDLT6   <- pDLT(proba6)
#pDLT6bis<-  pDLTbis(proba6)

# Output: Vector of mean(nTTB) and P(DLT) for the 6 doses
#--------------------------------------------------------
nTTB    <- c(mnTTB1,mnTTB2,mnTTB3,mnTTB4,mnTTB5,mnTTB6)
#print(paste("nTTB",sc,sep=""))
#paste(paste("nTTB",sep=""),sc) 
#print(paste("toto", 1:10,"=",x))
#paste("nTTB_sc",sc,sep="")
pDLT    <- c(pDLT1, pDLT2, pDLT3, pDLT4, pDLT5, pDLT6)
#pDLTbis <- c(pDLT1bis, pDLT2bis, pDLT3bis, pDLT4bis, pDLT5bis, pDLT6bis)
sc      <- rbind(listdo,nTTB,pDLT)



return(list(proba1=proba1,
            proba2=proba2,
            proba3=proba3,
            proba4=proba4,
            proba5=proba5,
            proba6=proba6,
            probaT1=probaT1,
            probaT2=probaT2,
            probaT3=probaT3,
            sc=sc)
      )
}
###Check hypothesis on the distribution of proba 
#sc1
#apply(desc_sc(1)$probaT1,2,which.max)
#apply(desc_sc(1)$probaT2,2,which.max)
#apply(desc_sc(1)$probaT3,2,which.max) 
#apply(desc_sc(1)$probaT1,1,sum)
#apply(desc_sc(1)$probaT2,1,sum)
#apply(desc_sc(1)$probaT3,1,sum)


#sc2
#apply(desc_sc(2)$probaT1,2,which.max)
#apply(desc_sc(2)$probaT2,2,which.max)
#apply(desc_sc(2)$probaT3,2,which.max)
#apply(desc_sc(2)$probaT1,1,sum)
#apply(desc_sc(2)$probaT2,1,sum)
#apply(desc_sc(2)$probaT3,1,sum)

#sc3  => pb
#apply(desc_sc(3)$probaT1,2,which.max)
#apply(desc_sc(3)$probaT2,2,which.max)
#apply(desc_sc(3)$probaT3,2,which.max)
#apply(desc_sc(3)$probaT1,1,sum)
#apply(desc_sc(3)$probaT2,1,sum)
#apply(desc_sc(3)$probaT3,1,sum)

#sc33  => ok
#apply(desc_sc(33)$probaT1,2,which.max)
#apply(desc_sc(33)$probaT2,2,which.max)
#apply(desc_sc(33)$probaT3,2,which.max)
#apply(desc_sc(33)$probaT1,1,sum)
#apply(desc_sc(33)$probaT2,1,sum)
#apply(desc_sc(33)$probaT3,1,sum)

#sc333  
#apply(desc_sc(333)$probaT1,2,which.max)
#apply(desc_sc(333)$probaT2,2,which.max)
#apply(desc_sc(333)$probaT3,2,which.max)
#apply(desc_sc(333)$probaT1,1,sum)
#round(apply(desc_sc(333)$probaT2,1,sum))
#apply(desc_sc(333)$probaT3,1,sum)

#sc4 =>ok 
#apply(desc_sc(4)$probaT1,2,which.max)
#apply(desc_sc(4)$probaT2,2,which.max)
#apply(desc_sc(4)$probaT3,2,which.max)
#apply(desc_sc(4)$probaT1,1,sum)
#apply(desc_sc(4)$probaT2,1,sum)
#apply(desc_sc(4)$probaT3,1,sum)

#sc5 
#apply(desc_sc(5)$probaT1,2,which.max)
#apply(desc_sc(5)$probaT2,2,which.max)
#apply(desc_sc(5)$probaT3,2,which.max)
#apply(desc_sc(5)$probaT1,1,sum)
#apply(desc_sc(5)$probaT2,1,sum)
#apply(desc_sc(5)$probaT3,1,sum)

#sc6 
#apply(desc_sc(6)$probaT1,2,which.max)
#apply(desc_sc(6)$probaT2,2,which.max)
#apply(desc_sc(6)$probaT3,2,which.max)
#apply(desc_sc(6)$probaT1,1,sum)
#apply(desc_sc(6)$probaT2,1,sum)
#apply(desc_sc(6)$probaT3,1,sum)

#sc7 
#apply(desc_sc(7)$probaT1,2,which.max)
#apply(desc_sc(7)$probaT2,2,which.max)
#apply(desc_sc(7)$probaT3,2,which.max)
#apply(desc_sc(7)$probaT1,1,sum)
#apply(desc_sc(7)$probaT2,1,sum)
#apply(desc_sc(7)$probaT3,1,sum)



#desc_sc(1)$sc
#write.table(desc_sc(1)$sc,"SC1.csv",sep=";")
#desc_sc(2)$sc
#write.table(desc_sc(2)$sc,"SC2.csv",sep=";")
#desc_sc(3)$sc
#write.table(desc_sc(5)$sc,"SC5.txt",sep="\t")
#write.table(desc_sc(6)$sc,"SC6.txt",sep="\t")
#write.table(desc_sc(7)$sc,"SC7.txt",sep="\t")

#print(xtable(desc_sc(4)$probaT2,digits=3),include.rownames=T,include.colomnames=T)
#print(xtable(desc_sc(4)$probaT3,digits=3),include.rownames=T,include.colomnames=T)

 #sc55 

#apply(desc_sc(55)$probaT1,2,which.max)
#apply(desc_sc(55)$probaT2,2,which.max)
#apply(desc_sc(55)$probaT3,2,which.max)


#desc_sc(55)$sc