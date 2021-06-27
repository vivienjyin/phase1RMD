 
\name{Tox2nTTP}
\alias{Tox2nTTP}
\title{Compute the nTTP score baed on the observed toxicity.}
\description{
  \code{Tox2nTTP} use the design parameters to compute the nTTP score baed on the observed toxicity.
}
\usage{

Tox2nTTP(tox = c(1,2,0),
 wm = matrix(c(0, 0.5, 0.75, 1 , 1.5,
 0, 0.5, 0.75, 1 , 1.5,
 0, 0 , 0 , 0.5, 1 ),
 byrow = T, ncol = 5),
 toxmax = 2.5)
}

\arguments{
  \item{tox}{Toxicity grades for the toxicity types.}
  \item{wm}{Clinical weight matrix, where toxicity types define the rows while the toxicity grades define the columns. Usually solicited from physicians.}
  \item{toxmax}{The normalization constant used in computing nTTP score. For details, see Ezzalfani et al(2013).}
}
\value{
  nTTP value.
}
\details{
The \code{Tox2nTTP} function calculates the nTTP values based on the weight matrix.

}

\examples{
#data(tox_matrix)

# Simulate the operating characteristics of clinical trials:
#Tox2nTTP(tox = c(1,2,0))
 
}
