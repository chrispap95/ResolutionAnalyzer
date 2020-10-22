#!/usr/bin/sh

# Basic configuration values - most are self explanatory
#     - filesToProcess: number of files to use for each dead fraction.
#     - splitLevel: maximum number of files per job.
#                   For no splitting, set it >= ${deadFractions[0]}
energyRange=0to3000
eta=1p7
pGenerator=SingleGamma
cmssw=${CMSSW_VERSION}
geometry=upgrade2023_D41
siteUrl=root://cmseos.fnal.gov/
deadFractions=(01 03 05 07)
filesToProcess=(298 99 60 43)
splitLevel=300
firstFile=0

# Checks if any dead fraction has been split
function isSplit() {
  for i in $@
  do
    if [ $i -ne 1 ]
    then
      echo 1
      exit
    fi
  done
  echo 0
}

# Prepare arrays for job splitting
configuration=${pGenerator}_E${energyRange}Eta${eta}
numberOfJobs=()
lastJobFiles=()
for files in ${filesToProcess[@]}
do
  if [ $((${files}%${splitLevel})) -eq 0 ]
  then
    numberOfJobs+=($((${files}/${splitLevel})))
    lastJobFiles+=(${splitLevel})
  else
    numberOfJobs+=($((${files}/${splitLevel} + 1)))
    lastJobFiles+=($((${files}%${splitLevel})))
  fi
done

# Print details about splitting
echo "Jobs are split as: (${numberOfJobs[@]})"
echo "Copy and paste this list to condorSubmitter_flat.sh"

# Loop over dead fractions and jobs per dead fraction
# nRuns: number of files to process
i=0
k=0
for df in ${deadFractions[@]}
do
  for j in `seq ${numberOfJobs[${k}]}`
  do
    if [ $(isSplit "${numberOfJobs[@]}") -eq 1 ]
    then
      namestring=E${energyRange}Eta${eta}_df${df}_${j}
    else
      namestring=E${energyRange}Eta${eta}_df${df}
    fi
    samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
    echo "outFilePath = out_${namestring}.root" > simpleBH_${namestring}.cfg
    echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> simpleBH_${namestring}.cfg
    echo "recoFileName = ntuples" >> simpleBH_${namestring}.cfg
    echo "deadfrac = 0.${df}" >> simpleBH_${namestring}.cfg
    echo "firstRun = ${firstFile}" >> simpleBH_${namestring}.cfg
    if [ ${j} -eq ${numberOfJobs[${k}]} ]
    then
      echo "nRuns = ${lastJobFiles[${k}]}" >> simpleBH_${namestring}.cfg
      firstFile=$((${firstFile} + ${lastJobFiles[${k}]}))
    else
      echo "nRuns = ${splitLevel}" >> simpleBH_${namestring}.cfg
      firstFile=$((${firstFile} + ${splitLevel}))
    fi
    ((i+=1))
  done
  ((k+=1))
done
