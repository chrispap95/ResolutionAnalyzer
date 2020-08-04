#!/bin/sh
USERBASE=`pwd`
rm ${CMSSW_VERSION}.tgz
cd ../../../
tar --exclude="*.root" --exclude=${CMSSW_BASE}/src/deadCellRegression --exclude-vcs -zcvf ${CMSSW_VERSION}.tgz ${CMSSW_VERSION}
mv ${CMSSW_VERSION}.tgz ${CMSSW_VERSION}/src/ResolutionAnalyzer
cd $USERBASE
