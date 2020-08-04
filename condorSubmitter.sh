#!/usr/bin/sh

for df in 01 03 05 07
do
for i in 5 10 15 20 30 40 60 80 100 140 200 280 400 550 750 1000 1400
do
cat > condor_E${i}Eta1p7_df${df}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_E${i}Eta1p7_df${df}.jdl
echo "Arguments = simpleBH_E${i}Eta1p7_df${df}.cfg out_E${i}Eta1p7_df${df}.root ${CMSSW_VERSION}" >> condor_E${i}Eta1p7_df${df}.jdl
cat >> condor_E${i}Eta1p7_df${df}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_E${i}Eta1p7_df${df}.jdl
done
done
