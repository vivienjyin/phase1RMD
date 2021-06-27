 
\name{plot.RunRMDVal}
\alias{plot.RunRMDVal}
\title{plot nTTP Toxicity-based in cycle 1 for all doses.}
\description{
  \code{plot.RunRMDVal} plot the dose-nTTP estimation in cycle based on the RunRMD estimation.
}

\usage{
\method{plot}{RunRMDVal}(x, \dots)
}

\arguments{
  \item{x}{The output of RunRMD. An S3 object of 'RunRMDVal' class.}
  \item{...}{Other plot parameters.}
}
\value{
  This function plots the dose-nTTP estimation in cycle based on the RunRMD estimation.
}
\details{
The \code{plot.RunRMDVal} function uses the patient data structure \code{patdata} to fit the linear mixed model as specified by \code{formula} based on the prior distributions as specified by \code{control}. The parameters of the MCMC simulation are specified by \code{iter}: the number of iterations, \code{burnin}: burn-ins, \code{thin}: thinning parameter and \code{chains}: number of chains in the MCMC simulation. The target nTTP score is defined by \code{tox.target}, default = 0.28.

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
 
# Load the example patient toxicity data:
#data(patdata)

# Estimate the dose-recommendation of the next patient cohort based on the toxicity data:
#res <- RunRMD(data=patdata, control=control, tox.target=0.28, 
#  iter=2000, burnin=500, thin=1, chains=1)

#plot(res);

}
