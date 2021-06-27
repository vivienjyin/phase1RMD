\name{GenToxProb}
\alias{GenToxProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Toxicity Probability Matrices Generation for Trial Simulation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to generate toxicity probability matrices employing proportional odds model, used as input for trial simulation. The generated matrices are a four-dimension array containing the probabilities of the occurrences of certain grades for certain types of toxicities, at each dose level and cycle under consideration. Dimension 1 refers to doses; dimension 2 corresponds to cycles of the treatment; dimension 3 regards the types of toxicities while dimenion 4 relates to grades.
}
\usage{
GenToxProb(toxtype = c("Neurological", "Renal", "Hematological"), 
           intercept.alpha = c(2, 3, 4.2, 5.7), 
           coef.beta = c(-0.2, -0.4, -0.7), 
           cycle.gamma = 0, MaxCycle = 6, Numdose = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toxtype}{A character vector, specifying toxicity types considered in the trial.}
  \item{intercept.alpha}{A four element numeric vector specifying the intercepts for the cumulative probabilities of the occurrences of grades 0-4 of toxicities in proportional odds model. See Details below.}
  \item{coef.beta}{A \eqn{n} numeric vector specifying the slope for dose in proportional odds model for \eqn{n} types of toxicities. See Details below.}
  \item{cycle.gamma}{A scalar controlling the cycle effect in simulation in proportional odds model. See Details below.}
  \item{MaxCycle}{Maximum of number of cycles of the treatment.}
  \item{Numdose}{The number of doses under investigation.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The probability matrices of the occurrences of certain grades for certain types of toxicities, at each dose level and cycle are simulated from a proportional odds model that takes \code{toxtype}, \code{intercept.alpha}, \code{coef.beta} and \code{cycle.gamma} as the inputs. The model is as follows:

\deqn{logit(c_j) = \alpha_j + \beta * dose + \gamma * cycle, j=0,1,2,3}
where given a dose and cycle, \eqn{c_j} is the cumulative probabilities of occurrences of toxicities grades up to \eqn{j}; \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3} are the intercepts, making up \code{intercept.alpha}; \code{coef.beta} stores \eqn{\beta} for different types of toxicities while \code{cycle.gamma} acts as \eqn{\gamma}, controlling the cycle effect. The probability for the occurrence of a certain grade of a certain toxicity type at a given dose and cycle is obtained from taking the corresponding differences between cumulative probabilities.
}
\value{
A four-dimension array where dimension 1 refers to doses; dimension 2 corresponds to cycles of the treatment; dimension 3 regards the types of toxicities while dimenion 4 relates to grades.
}

\examples{
#tox.matrix <- GenToxProb()
}
