 
\name{SimRMD}
\alias{SimRMD}
\title{Simulate trial characters of the phase I dose-finding design.}
\description{
  SimRMD uses the design parameters to simulate the trial characters of phase I dose-finding design.
}
\usage{
SimRMD(seed=2014, strDose=1, chSize=3, trlSize=36,
   numTrials=1000, sdose=1:6, MaxCycle=6, tox.target=0.28,
   control, iter=10000, burnin=4000, thin=1, chains=1,    
   pathout='./', tox.matrix, wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0.5, 0.75, 1  , 1.5, 
                                        0, 0  , 0   , 0.5, 1  ), 
                                        byrow = T, ncol = 5), toxmax = 2.5)
}

\arguments{
  \item{seed}{Random seed of the simulation (default = 2014).}
  \item{strDose}{Start dose (default = 1).}
  \item{chSize}{Size of each patient cohort (default = 3).}
  \item{trlSize}{Total number of patients in the trial (default = 36).}
  \item{numTrials}{Number of simulated trials (default = 1000).}
  \item{sdose}{Dose range (default = 1:6).}
  \item{MaxCycle}{Maximum treatment cycle (default = 6).}
  \item{tox.target}{Target toxicity score (default = 0.28).}
  \item{control}{Prior distributions as specified by \code{parm} function.}
  \item{iter}{Total number of MCMC simulations (default = 4000).}
  \item{burnin}{Number of burn-ins in the MCMC simulation (default = 1000).}
  \item{thin}{Thinning parameter (default = 1).}
  \item{chains}{Number of chains in the MCMC simulation (default = 1).} 
  \item{pathout}{The output directory of simulation results (default = './').}
  \item{tox.matrix}{Four-dimension array contains Dose, Cycle, Toxicity Type, Toxicity Grade.}
  \item{wm}{Clinical weight matrix, where toxicity types define the rows while the toxicity grades define the columns. Usually solicited from physicians.}
  \item{toxmax}{The normalization constant used in computing nTTP score. For details, see Ezzalfani et al(2013).}
}

\value{
  A list containing the following elements:
  \item{DOSE-RECOMMENDED}{Recommended dose for the next patient cohort. }
  \item{Estimate}{Mean, SD, and median estimates of the toxicity score at each dose level.}
  \item{Quantiles}{2.5\%, 25\%, 50\%, 75\%, 97.5\% quantile estimates of the toxicity score at each dose level.}
}

\details{
The \code{SimRMD} function uses the patient data structure \code{patdata} to fit the linear mixed model as specified by \code{formula} based on the prior distributions as specified by \code{control}. The parameters of the MCMC simulation are specified by \code{iter}: the number of iterations, \code{burnin}: burn-ins, \code{thin}: thinning parameter and \code{chains}: number of chains in the MCMC simulation. The target nTTP score is defined by \code{tox.target}, default = 0.28.
}

\examples{
# Setup the prior distributions for the phase I dose-finding model:
control <- list(
 beta.dose = parm("normal", mean = 0, var = 1000),
 beta.other = parm("normal", mean = 0, var = 1000 ),
 gamma = parm("normal", mean = 0, var = 100 ),
 s2.gamma = parm("invgamma", shape = 0.001, scale = 0.001),
 s2.epsilon = parm("invgamma", shape = 0.001, scale = 0.001)
 )


# Generate the toxicity matrix
# tox.matrix <- GenToxProb(
#    toxtype = c("Renal", "Neuro", "Heme"), 
#	intercept.alpha = c(2, 3, 4.2, 5.7), 
#	coef.beta = c(-0.2, -0.4, -0.7), 
#	cycle.gamma = 0)

# Simulate the dose-recommendation based on the toxicity matrix
# simu <- SimRMD(seed=2014, strDose=1, chSize=3, trlSize=12,
#    numTrials=1, sdose=1:6, MaxCycle=5, tox.target=0.28,
#    control=control, iter=10, burnin=2, thin=1, chains=1,    
#    pathout='./', tox.matrix=tox.matrix)
	
}
