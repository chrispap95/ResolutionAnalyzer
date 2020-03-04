#!/bin/sh
USERBASE=`pwd`
rm CMSSW_10_6_3_patch1.tgz
cd ../../../
tar --exclude="*.root" --exclude=${CMSSW_BASE}/src/deadCellRegression --exclude-vcs -zcvf CMSSW_10_6_3_patch1.tgz CMSSW_10_6_3_patch1/
# xrdcp -f CMSSW_10_6_3_patch1.tgz root://cmseos.fnal.gov//store/user/chpapage/CMSSW_10_6_3_patch1.tgz
mv CMSSW_10_6_3_patch1.tgz CMSSW_10_6_3_patch1/src/ResolutionAnalyzer
cd $USERBASE
