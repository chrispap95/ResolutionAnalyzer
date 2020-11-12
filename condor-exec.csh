#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
tar -xf ${3}.tgz
rm ${3}.tgz
setenv SCRAM_ARCH slc7_amd64_gcc820
cd ${3}/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd ResolutionAnalyzer
mkdir {lib,bin,obj}
make
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.67.0/lib/
#echo "Arguments passed to this script are: for 1: $1, and for 2: $2"
./bin/simpleBH -c scripts/$1

# Export file
set nonomatch
if ( $1 == *"to"* ) then
xrdcp -f $2 root://cmseos.fnal.gov//store/user/${4}/DeadCellsSamples/TrainingSamples/$2
else
xrdcp -f $2 root://cmseos.fnal.gov//store/user/${4}/DeadCellsSamples/EvaluationSamples/$2
endif
### remove the output file if you don't want it automatically transferred when the job ends
rm $2
cd ${_CONDOR_SCRATCH_DIR}
rm -rf $3
