#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
### for case 1. EOS have the following line, otherwise remove this line in case 2.
#xrdcp -s root://cmseos.fnal.gov//store/user/chpapage/CMSSW_10_6_3_patch1.tgz .
tar -xf CMSSW_10_6_3_patch1.tgz
rm CMSSW_10_6_3_patch1.tgz
setenv SCRAM_ARCH slc7_amd64_gcc820
cd CMSSW_10_6_3_patch1/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd ResolutionAnalyzer
mkdir {lib,bin,obj}
make
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:lib
#echo "Arguments passed to this script are: for 1: $1, and for 2: $2"
./bin/simpleBH -c scripts/$1
xrdcp -f $2 root://cmseos.fnal.gov//store/user/chpapage/$2
### remove the output file if you don't want it automatically transferred when the job ends
rm $2
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_6_3_patch1
