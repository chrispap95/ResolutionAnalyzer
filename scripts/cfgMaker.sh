#!/usr/bin/sh

# Basic configuration values - self explanatory

eta=1p7
pGenerator=SingleGamma
cmssw=CMSSW_10_6_3_patch1
geometry=upgrade2023_D41

# Loop over dead fractions and energies
# nRuns: number of files to process

for df in 01 03 05 07
do
for En in 10 20 50 100
do
namestring=E${En}Eta${eta}_df${df}
configuration=${pGenerator}_E${En}Eta${eta}
samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
echo "outFilePath = out_${namestring}.root" > simpleBH_${namestring}.cfg
echo "filePath = root://cmseos.fnal.gov//${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> simpleBH_${namestring}.cfg
cat >> simpleBH_${namestring}.cfg << "EOF"
recoFileName = ntuples
nRuns = 200
EOF
echo "deadfrac = 0.${df}" >> simpleBH_${namestring}.cfg
done
done
