# phase1RMD
phase1RMD R package

```
> #load the patient efficacy data
> data(eff_dat)
> #load the patient toxicity data
> data(tox_dat)
```

```
> RunRMDEFF(efficacy.dat = eff_dat, toxicity.dat = tox_dat)
Model : RMD with longitudinal toxicity


Doses(skeleton):
 	1 	2 	3 	4 	5 	6 


We are recommending the dose for your next cohort of patients...
The maximum sample size is : 36
The current enrolled number of patients are: 33
The current enrolled cohort is: 11
You are right now in stage 2, randomizing the next cohort of patients towards higher predicted efficacy...
Posterior estimates (mean) of toxicity and efficacy:
                       dose1      dose2     dose3     dose4     dose5     dose6
toxicity_profile1 0.02294106 0.08046547 0.1379899 0.1955143 0.2530387 0.3105631
toxicity_profile2 0.02244048 0.07996489 0.1374893 0.1950137 0.2525381 0.3100625
efficacy_profile  0.05793578 0.18123673 0.3351263 0.5196046 0.7346715 0.9803271
Next recommended dose: 5
$nxtdose
[1] 5

$tox.pf
                       dose1      dose2     dose3     dose4     dose5     dose6
toxicity_profile1 0.02294106 0.08046547 0.1379899 0.1955143 0.2530387 0.3105631
toxicity_profile2 0.02244048 0.07996489 0.1374893 0.1950137 0.2525381 0.3100625

$eff.pf
                      dose1     dose2     dose3     dose4     dose5     dose6
efficacy_profile 0.05793578 0.1812367 0.3351263 0.5196046 0.7346715 0.9803271

$allow.doses
[1] 1 2 3 4 5


$eff.pf
                      dose1     dose2     dose3     dose4     dose5     dose6
efficacy_profile 0.05793578 0.1812367 0.3351263 0.5196046 0.7346715 0.9803271

$allow.doses
[1] 1 2 3 4 5
```
