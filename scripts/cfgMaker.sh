#!/usr/bin/sh

# Basic configuration values - self explanatory

eta=1p7
pGenerator=SingleGamma
cmssw=CMSSW_10_6_3_patch1
geometry=upgrade2023_D41
siteUrl=root://cmseos.fnal.gov/
deadFractions=(01 03 05 07)
energies=(5 10 20 40 60 80 100)
filesToProcess=200

# Loop over dead fractions and energies
# nRuns: number of files to process

for df in ${deadFractions[@]}
do
for En in ${energies[@]}
do
namestring=E${En}Eta${eta}_df${df}
configuration=${pGenerator}_E${En}Eta${eta}
samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
echo "outFilePath = out_${namestring}.root" > simpleBH_${namestring}.cfg
echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> simpleBH_${namestring}.cfg
echo "recoFileName = ntuples" >> simpleBH_${namestring}.cfg
echo "nRuns = ${filesToProcess}" >> simpleBH_${namestring}.cfg
echo "deadfrac = 0.${df}" >> simpleBH_${namestring}.cfg
done
done
