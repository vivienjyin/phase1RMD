\name{RunRMDEFF}
\alias{RunRMDEFF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive Implementation for an Adaptive, Multi-Stage Phase I Dose-Finding Design 
}
\description{
A function to interactively implement an adaptive, multi-stage phase 1 dose-finding design incorporating continuous efficacy and toxicity data from multiple treatment cycles, proposed by Du et al(2017)
}
\usage{
RunRMDEFF(efficacy.dat = NULL, toxicity.dat, trialSize = 36,
          seed = 624, chSize = 3, MaxCycle = 6, doses = 1:6, 
          tox.target = 0.28, p1 = 0.2, p2 = 0.2, ps1 = 0.2, 
          thrd1 = 0.28, thrd2 = 0.28, proxy.thrd = 0.1, 
          dose_flag = 0, param.ctrl = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{efficacy.dat}{A data frame specifying the efficacy data, which must be of the following form. The first column is the subject ID (character), named "subID", where "subID" consists of a cohort number and a subject number. For example, "cohort2subject3" denotes the third subject in the second cohort; the second column is dose allocated (integer), named "dose" and the third column is the continuous efficacy (numeric), named "Efficacy". This data frame is sorted by "subID". This argument is set to NULL if none of the subjects efficacy measurement is available.}
  \item{toxicity.dat}{A data frame specifying the toxicity data from multiple treatment cycles, which must be of the following form. The first column is the subject ID (character), named "subID", where "subID" consists of a cohort number and a subject number. For example, "cohort2subject3" denotes the third subject in the second cohort; the second column is dose allocated (integer), named "dose"; the third column, "cycle", denotes the cycle information (integer); "nTTP", the fourth column shows the corresponding nTTP score and the fifth column is binary DLT event indicator (integer), named "DLT". This data frame is sorted first by cohort number, second by "cycle" and then by subject number.}
  \item{trialSize}{The maximum sample size for the trial under investigation. Default is 36. Must be the multiple of cohort size, represented by \code{chSize}.}
  \item{seed}{A number to make simulation results reproducible.}
  \item{chSize}{The cohort size of patients recruited. Default is 3.}
  \item{MaxCycle}{The maximum number of cycles for the treatment, Default is 6.}
  \item{doses}{A vector of doses that users are going to explore. Default is \code{1:6}, where dose 1 through dose 6 are being tested.}
  \item{tox.target}{The target toxicity at cycle 1 of the treatment. Default is 0.28.}
  \item{p1}{The probability cutoff for cycle 1 toxicity. Default is 0.2. See Details below.}
  \item{p2}{The probability cutoff for late cycles toxicity beyond cycle 1. Default is 0.2. See Details below.}
  \item{ps1}{The probability cutoff for defining allowable (safe) doses in stage 1. Default is 0.2. See Details below.}
  \item{thrd1}{An upper bound of toxicity for cycle 1 of the treatment. Default is 0.28. See Details below.}
  \item{thrd2}{An upper bound of toxicity for late cycles of the treatment, beyond cycle 1. Default is 0.28. See Details below.}
  \item{proxy.thrd}{A distance parameter used to define efficacious doses. Any dose whose predicted efficacy is within \code{proxy.thrd} away from the largest one among the safe doses will be declared an efficacious dose.}
  \item{dose_flag}{A binary indicator. Value 1 indicates that dose has not escalated to dose 2 since the first enrolled cohort while value 0 shows otherwise. Note that for the first enrolled cohort, it is set at 0. The first cohort of patients is treated at the starting dose, dose 1.}
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
\code{p1},\code{p2},\code{thrd1} and \code{thrd2} are used to define allowable (safe) doses through the two probability conditions: \deqn{P(tox1 < thrd1) > p1}{P(tox1 < thrd1) > p1}\deqn{P(tox2 < thrd2) > p2}{P(tox2 < thrd2) > p2}, where tox1 and tox2 denote the estimated toxicity at cycle 1 of the treament and late cycles, respectively. In stage 1, \code{ps1} acts as \code{p1} and \code{p2}.

Note that our design allows patients dropout due to DLT events, so the data frames \code{efficacy.dat} and \code{toxicity.dat} only include the information for the available patients and available cycles.
}



\value{
The function prints information regarding the posterior toxicity/efficacy profile for each dose under exploration, the recommended dose for the next cohort of patients, and return a list with the following components:
\item{nxtdose}{The next recommended dose. For dose recommendation at the end of the trial, see \code{opt.dose}.}
\item{tox.pf}{A data frame demonstrating the toxicity profile for each dose under exploration, for cycle 1 of the treatment (row 1) and the other late cycles beyond cycle 1 as an averaged measure (row 2).}
\item{eff.pf}{A data frame demonstrating the efficacy profile for each dose under exploration, reported after stage 2 of the design begins.}
\item{allow.doses}{A vector of allowable (safe) doses, defined through the two probability conditions, specified in Details.}
\item{dose_flag}{A binary indicator. Value 1 indicates that dose has not escalated to dose 2 since the first enrolled cohort while value 0 shows otherwise. Used as input for the next cohort of patients.}
\item{opt.dose}{The dose that is recommended at the end of the trial. We recommend the lowest (safest) dose that is efficacious. The physicians, however,  may make their own desicions by reviewing the efficacy and toxicitiy profile, returned by \code{eff.pf} and \code{tox.pf}, respectively.}
\item{p_nTTP1}{Boxplots showing posterior nTTP for cycle 1 of the treatment across dose levels under investigation.}
\item{p_nTTP2}{Boxplots showing posterior nTTP for late cycles of the treatment across dose levels under investigation.}
\item{p_EFF}{Boxplots showing posterior efficacy of the treatment across dose levels under investigation, reported after stage 2 of the design begins.}
}

\examples{

#load the patient efficacy data
#data(eff_dat)
#head(eff_dat)
#----------------------#
#load the patient toxicity data
#data(tox_dat)
#head(tox_dat)

#RunRMDEFF(efficacy.dat = eff_dat, toxicity.dat = tox_dat)

}
