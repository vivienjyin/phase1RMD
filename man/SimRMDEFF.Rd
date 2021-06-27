\name{SimRMDEFF}
\alias{SimRMDEFF}
\title{
Simulation for an Adaptive, Multi-Stage Phase I Dose-Finding Design 
}
\description{
A function to implement simulation for an adaptive, multi-stage phase 1 dose-finding design incorporating a continuous efficacy outcome and toxicity data from multiple treatment cycles, proposed by Du et al(2017)
}
\usage{
SimRMDEFF(numTrials = 100, trialSize = 36, doses = 1:6, cycles = 1:6,  
          eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9), 
          eff.sd = 0.2, tox.target = 0.28, p1 = 0.2, p2 = 0.2, 
          ps1 = 0.2, StrDose = 1, chSize = 3, 
          tox.matrix = NULL, 
          proxy.thrd = 0.1, thrd1 = 0.28, thrd2 = 0.28, 
          wm = matrix(c(0, 0.5, 0.75, 1  , 1.5, 
                        0, 0.5, 0.75, 1 , 1.5, 
                        0, 0  , 0   , 0.5, 1  ), 
                        byrow = T, ncol = 5), 
          toxmax = 2.5, toxtype = NULL, intercept.alpha = NULL, 
          coef.beta = NULL, cycle.gamma = NULL, 
          param.ctrl = list())
}
\arguments{
  \item{numTrials}{An integer specifying the number of simulations.}
  \item{trialSize}{The maximum sample size for trial simulation. Default is 36. Must be the multiple of cohort size, represented by \code{chSize}.}
  \item{doses}{A vector of doses that users are going to explore. Default is \code{1:6}, where dose 1 through dose 6 are being tested.}
  \item{cycles}{A vector of cycles that the treatment plans to go through. Default is \code{1:6}, where patients will experience up to 6 cycles of the treatment.}
  \item{eff.structure}{A vector of the mean values of efficacy outcome for each dose level explored. Default is \code{c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9)}, corresponding to dose 1 through dose 6, where the efficacy increases with dose.}
  \item{eff.sd}{The standard deviation for generating the efficacy outcome for each patient. Default is 0.2, which already gives large skewness/variation.}
  \item{tox.target}{The target toxicity at cycle 1 of the treatment. Default is 0.28.}
  \item{p1}{The probability cutoff for cycle 1 toxicity. Default is 0.2. See Details below.}
  \item{p2}{The probability cutoff for later cycles toxicity beyond cycle 1. Default is 0.2. See Details below.}
  \item{ps1}{The probability cutoff for defining allowable (safe) doses in stage 1. Default is 0.2. See Details below.}
  \item{StrDose}{A number to specify the starting dose. Default is 1.}
  \item{chSize}{The cohort size of patients recruited. Default is 3.}
  \item{tox.matrix}{Optional. A four-dimension array specifying the probabilities of the occurrences of certain grades for certain types of toxicities, at each dose level and cycle under consideration. Dimension 1 refers to doses; dimension 2 corresponds to cycles of the treatment; dimension 3 regards the types of toxicities while dimenion 4 relates to grades. If null, which is default choice, the arguments \code{toxtype}, \code{intercept.alpha}, \code{coef.beta}, \code{cycle.gamma} must be provided to simulate this array.}
  \item{proxy.thrd}{A distance parameter to define efficacious doses. Any dose whose predicted efficacy is within \code{proxy.thrd} away from the largest one among the safe doses will be declared an efficacious dose.}
  \item{thrd1}{An upper bound of toxicity for cycle 1 of the treatment. Default is 0.28. See Details below.}
  \item{thrd2}{An upper bound of toxicity for late cycles of the treatment, beyond cycle 1. Default is 0.28. See Details below.}
  \item{wm}{Clinical weight matrix, where toxicity types define the rows while the toxicity grades define the columns. Usually solicited from physicians.}
  \item{toxmax}{The normalization constant used in computing nTTP score. For details, see Ezzalfani et al(2013).}
  \item{toxtype}{Only specified when \code{tox.matrix} is null. This argument, a character vector, specifies toxicity types considered in the trial.}
  \item{intercept.alpha}{Only specified when \code{tox.matrix} is null. A four element numeric vector specifying the intercepts for the cumulative probabilities of the occurrences of grades 0-4 of toxicities in proportional odds model. See Details below.}
  \item{coef.beta}{Only specified when \code{tox.matrix} is null. A \eqn{n} numeric vector specifying the slope for dose in proportional odds model for \eqn{n} types of toxicities. See Details below.}
  \item{cycle.gamma}{Only specified when \code{tox.matrix} is null. A scalar controlling the cycle effect in simulation in proportional odds model. See Details below.}
  \item{param.ctrl}{A list specifying the prior distribution for the parameters. \code{p1_beta_intercept}, the prior mean of intercept of toxicity model assuming a normal prior;
\code{p2_beta_intercept}, the precision (inverse of variance) of intercept of toxicity model assuming a normal prior;
\code{p1_beta_cycle}, the prior mean of cycle effect of toxicity model assuming a normal prior; 
\code{p2_beta_cycle}, the precision (inverse of variance) of cycle effect of toxicity model assuming a normal prior;
\code{p1_beta_dose}, the prior minimum of dose effect of toxicity model assuming a uniform prior;  
\code{p2_beta_dose}, the prior maximum of dose effect of toxicity model assuming a uniform prior;  
\code{p1_alpha}, the prior mean vector of the parameters from efficacy model assuming a multivariate normal prior;
\code{p2_alpha}, the prior precision matrix (inverse of covariance matrix) of the parameters from efficacy model assuming a multivariate normal prior;
\code{p1_gamma0}, the prior mean of association parameter \eqn{\gamma} (See Du et al(2017)) of two submodels of the joint model assuming a normal prior;
\code{p2_gamma0}, the prior precision (inverse of variance) of association parameter \eqn{\gamma} of two submodels of the joint model assuming a normal prior. Default is non-informative priors.}
}

\details{
\code{p1},\code{p2},\code{thrd1} and \code{thrd2} are used to define allowable (safe) doses through the two probability conditions: \deqn{P(tox1 < thrd1) > p1}{P(tox1 < thrd1) > p1}\deqn{P(tox2 < thrd2) > p2}{P(tox2 < thrd2) > p2}, where tox1 and tox2 denote the toxicity at cycle 1 of the treament and later cycles, respectively. In stage 1, \code{ps1} acts as \code{p1} and \code{p2}.

The user can explore different dose-efficacy patterns using the argument \code{eff.structure}. For example, \code{eff.structure = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)} corresponds to a flat relationship where the efficacy does not change with dose.

When \code{tox.matrix} is null, the probability matrices of the occurrences of certain grades for certain types of toxicities, at each dose level and cycle will be simulated from a proportional odds model that takes \code{toxtype}, \code{intercept.alpha}, \code{coef.beta} and \code{cycle.gamma} as the inputs. The model is as follows:

\deqn{logit(c_j) = \alpha_j + \beta * dose + \gamma * cycle, j=0,1,2,3}
where given a dose and cycle, \eqn{c_j} is the cumulative probabilities of occurrences of toxicities grades up to \eqn{j}; \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3} are the intercepts, making up \code{intercept.alpha}; \code{coef.beta} stores \eqn{\beta} for different types of toxicities while \code{cycle.gamma} acts as \eqn{\gamma}, controlling the cycle effect. The probability for the occurrence of a certain grade of a certain toxicity type at a given dose and cycle is obtained from taking the corresponding differences between cumulative probabilities.

}


\value{
A list with the following components:
\item{op.table}{A table reporting the operating characteristics of the trial simulation.}
\item{sc}{The assumed true mean longitudinal toxicity and efficacy structure for each dose level under investigation, along with the probability of DLT event, for cycle 1 of the treatment.}
}


\examples{

#######Dose 5 target toxicity dose, flat cycle effect#####################

###############a flat dose-efficacy pattern###############################

# simul <- SimRMDEFF(toxtype = c("H", "L", "M"), 
                   # intercept.alpha = c(1.9, 2.3, 2.6, 3.1), 
                   # coef.beta = c(-0.3, -0.2, -0.25), 
                   # cycle.gamma = 0, tox.target = 0.23, 
                   # thrd1 = 0.23, thrd2 = 0.23,
                   # eff.structure = rep(0.5, 6),
                   # p1 = 0.1, p2 = 0.1, ps1 = 0.1, 
                   # numTrials = 1)


}
