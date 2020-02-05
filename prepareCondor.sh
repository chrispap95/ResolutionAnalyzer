#!/bin/sh
USERBASE=`pwd`
cd ../../../
tar --exclude="*.root" --exclude-vcs -zcvf CMSSW_10_6_4.tgz CMSSW_10_6_4/
# xrdcp -f CMSSW_10_6_3_patch1.tgz root://cmseos.fnal.gov//store/user/chpapage/CMSSW_10_6_3_patch1.tgz
mv CMSSW_10_6_4.tgz CMSSW_10_6_4/src/ResolutionAnalyzer
cd $USERBASE
