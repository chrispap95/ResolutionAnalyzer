# ResolutionAnalyzer
Code that performs simpleBH analysis using CMSSW produced ntuples as input.

## Introduction
Works on ```CMSSW_10_6_*```. Setup an environment and clone into ```$CMSSW_BASE/src```.
Then, do
```bash
mkdir {bin,lib,obj}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`pwd`/lib
make
```

## Execute code
To execute the code
```bash
./bin/simpleBH -c script/simpleBH.cfg
```

You can view the output
```bash
root -l out.root
```
