#!/usr/bin/sh
source prepareCondor.sh

eta=1p7

for df in 01 03 05 07
do
for En in 5 10 15 20 30 40 60 80 100 140 200 280 400 550 750 1000 1400
do
namestring=E${En}Eta${eta}_df${df}
cat > condor_${namestring}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_${namestring}.jdl
echo "Arguments = simpleBH_${namestring}.cfg out_${namestring}.root ${CMSSW_VERSION}" >> condor_${namestring}.jdl
cat >> condor_${namestring}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${namestring}.jdl
done
done
