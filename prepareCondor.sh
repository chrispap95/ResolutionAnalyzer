#!/bin/sh
USERBASE=`pwd`
rm ${CMSSW_VERSION}.tgz
cd ../../../
echo "Creating tarball..."
tar --exclude="*.root" --exclude=${CMSSW_BASE}/src/deadCellRegression --exclude-vcs -zcf ${CMSSW_VERSION}.tgz ${CMSSW_VERSION}
mv ${CMSSW_VERSION}.tgz ${CMSSW_VERSION}/src/ResolutionAnalyzer
cd $USERBASE
if [ ! -f ${CMSSW_VERSION}.tgz ]; then
echo "Error: tarball doesn't exist!"
else
echo -n " Done!"
fi 
