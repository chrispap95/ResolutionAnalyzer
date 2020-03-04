#!/usr/bin/sh

for i in `seq 0 19`
do
for E in 400 550 750 1000 1400
do
for df in 7
do
cat > condor_E${E}Eta1p7_df0${df}_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = condor-exec.csh, CMSSW_10_6_3_patch1.tgz
EOF
echo "Arguments = simpleBH_E${E}Eta1p7_df0${df}_${i}.cfg out_E${E}Eta1p7_df0${df}_${i}.root" >> condor_E${E}Eta1p7_df0${df}_${i}.jdl
cat >> condor_E${E}Eta1p7_df0${df}_${i}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_E${E}Eta1p7_df0${df}_${i}.jdl
done
done
done
