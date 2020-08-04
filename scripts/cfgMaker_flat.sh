#!/usr/bin/sh

# Basic configuration values - self explanatory

eta=1p7
pGenerator=SingleGamma
cmssw=CMSSW_10_6_3_patch1
geometry=upgrade2023_D41
namestring=E0to3000Eta${eta}_df01
configuration=${pGenerator}_E0to3000Eta${eta}
samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/

# Loop over dead fractions and energies
# nRuns: number of files to process

for i in `seq 0 499`
do
echo "outFilePath = out_${namestring}_${i}.root" > simpleBH_${namestring}_${i}.cfg
echo "filePath = root://cmseos.fnal.gov//${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> simpleBH_${namestring}_${i}.cfg
cat >> simpleBH_${namestring}_${i}.cfg << "EOF"
recoFileName = ntuples
nRuns = 1
EOF
echo "firstRun = $(($i * 1))" >> simpleBH_${namestring}_${i}.cfg
echo "deadfrac = 0.01" >> simpleBH_${namestring}_${i}.cfg
done
