 
\name{parm}
\alias{parm}
\title{ Create prior for the phase 1 dose-finding model}
\description{
  \code{parm} creates an S3 object of class \code{parm} that contains all necessary prior information to estimate the phase I repeated measures model.
}
\usage{
parm(prior = c("gamma", "invgamma", "normal"), 
   mean=0, var=100, shape=0.001, scale=0.001)
}
\arguments{
  \item{prior}{ a symbolic description of the prior distribution (gamma, invgamma, normal). }
  \item{mean}{ the prior parameters that are dependent on the distribution that is used. Specifications are gamma(shape, rate), invgamma(shape, scale) and normal(mean,var). }
  \item{var}{ the prior parameters that are dependent on the distribution that is used. Specifications are gamma(shape, rate), invgamma(shape, scale) and normal(mean,var). }
  \item{shape}{ the prior parameters that are dependent on the distribution that is used. Specifications are gamma(shape, rate), invgamma(shape, scale) and normal(mean,var). }
  \item{scale}{ the prior parameters that are dependent on the distribution that is used. Specifications are gamma(shape, rate), invgamma(shape, scale) and normal(mean,var). }
}
\value{
  An object of class \code{parm}, containing information for the prior.
}
\details{
  The function \code{parm} builds three types of priors: gamma, invgamma and normal. 

  Gamma parameters:
  \itemize{
    \item \bold{shape} : shape parameter for Gamma prior for sigma
    \item \bold{rate} : rate parameter for Gamma prior for sigma
  }

  Inverse gamma parameters:
  \itemize{
    \item \bold{shape} : shape parameter for the inverse gamma prior on sigma.
    \item \bold{scale} : scale parameter for the inverse gamma prior on sigma.
  }

  Normal parameters:
  \itemize{
    \item \bold{mean} : mean parameter for the normal prior.)
    \item \bold{var} : variance parameter for the normal prior.)
  }

}
\author{
  Jun Yin
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
}
