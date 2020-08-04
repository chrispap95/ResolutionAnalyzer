#!/usr/bin/sh

for i in `seq 0 0`
do
cat > condor_E0to3000Eta1p7_df01_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_E0to3000Eta1p7_df01_${i}.jdl
echo "Arguments = simpleBH_E0to3000Eta1p7_df01_${i}.cfg out_E0to3000Eta1p7_df01_${i}.root" >> condor_E0to3000Eta1p7_df01_${i}.jdl
cat >> condor_E0to3000Eta1p7_df01_${i}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_E0to3000Eta1p7_df01_${i}.jdl
done
