#!/usr/bin/sh

# Basic configuration values - self explanatory

energyRange=0to3000
eta=1p7
pGenerator=SingleGamma
cmssw=CMSSW_10_6_3_patch1
geometry=upgrade2023_D41
siteUrl=root://cmseos.fnal.gov/
deadFractions=(01 03 05 07)
filesToProcess=(298 99 60 43)
firstFile=0
configuration=${pGenerator}_E${energyRange}Eta${eta}

# Loop over dead fractions and energies
# nRuns: number of files to process

i=0
for df in ${deadFractions[@]}
do
namestring=E${energyRange}Eta${eta}_df${df}
samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
echo "outFilePath = out_${namestring}_${i}.root" > simpleBH_${namestring}_${i}.cfg
echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> simpleBH_${namestring}.cfg
echo "recoFileName = ntuples" >> simpleBH_${namestring}.cfg
echo "nRuns = ${filesToProcess[$i]}" >> simpleBH_${namestring}.cfg
echo "deadfrac = 0.${df}" >> simpleBH_${namestring}.cfg
echo "firstRun = ${firstFile}" >> simpleBH_${namestring}_${i}.cfg
firstFile=$(($firstFile + $filesToProcess[$i]))
((i+=1))
done
