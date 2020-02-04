# ResolutionAnalyzer
Code that performs simpleBH analysis using CMSSW produced ntuples as input.

## Introduction
Works on ```CMSSW_10_6_*```. Setup an environment and clone into ```$CMSSW_BASE/src```.
Then, do
```bash
mkdir {bin,lib,obj}
source setup.sh
make
```

## Execute code
To execute the code
```bash
./bin/simpleBH -c scripts/simpleBH.cfg
```

You can view the output
```bash
root -l out.root
```

## Use condor to submit
*__Beware__: These instructions and scripts currently work only on the LPC cluster and contain my username. In order to use them, you have to put your username and if you are not on the LPC cluster you need to edit the code accordingly.*

I strongly suggest that you setup a fresh CMSSW release in order to make Condor submissions. So, do
```bash
mkdir condorSubmissions
cd condorSubmissions
cmsrel CMSSW_10_6_3_patch1
cd CMSSW_10_6_3_patch1/src
cmsenv
git clone https://github.com/chrispap95/ResolutionAnalyzer.git
scram b
```
Then, prepare your CMSSW to trasfer it to the nodes by issuing
```bash
cd ResolutionAnalyzer
sh condorPrepare.sh
```

Then, submit a job
```bash
condor_submit condor.jdl
```
