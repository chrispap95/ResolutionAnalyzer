#!/usr/bin/sh
source prepareCondor.sh

eta=1p7
namestring=E0to3000Eta${eta}_df01

for i in `seq 0 0`
do
cat > condor_${namestring}_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_${namestring}_${i}.jdl
echo "Arguments = simpleBH_${namestring}_${i}.cfg out_${namestring}_${i}.root" >> condor_${namestring}_${i}.jdl
cat >> condor_${namestring}_${i}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${namestring}_${i}.jdl
done
