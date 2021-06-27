# phase1RMD
phase1RMD R package
Author: Jun (Vivien) Yin

## Load an example of the patient efficacy data
```
> #load the patient efficacy data
> data(eff_dat)
> head(eff_dat)
            subID dose    Efficacy
1 cohort1subject1    1 0.014692986
2 cohort1subject2    1 0.005370450
3 cohort1subject3    1 0.004324666
4 cohort2subject1    2 0.005531986
5 cohort2subject2    2 0.300249297
6 cohort2subject3    2 0.002631852
```

## Load an example of the patient toxicity data
```
> #load the patient toxicity data
> data(tox_dat)
> head(tox_dat)
            subID dose cycle nTTP DLT
1 cohort1subject1    1     1  0.0   0
2 cohort1subject2    1     1  0.0   0
3 cohort1subject3    1     1  0.0   0
4 cohort1subject1    1     2  0.2   0
5 cohort1subject2    1     2  0.0   0
6 cohort1subject3    1     2  0.0   0
```

## Load an example of the patient toxicity data
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
