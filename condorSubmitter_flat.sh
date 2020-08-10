#!/usr/bin/sh
source prepareCondor.sh

# Define submission parameters
#   - samplesNumber is the number of files to process per dead fraction
#     to find the proper number you need to weight the samples with the
#     dead fractions.
energyRange=0to3000
eta=1p7
deadFractions=(01 03 05 07)

for df in ${deadFractions[@]}
do
namestring=E${energyRange}Eta${eta}_df${df}
argument=simpleBH_${namestring}.cfg\ out_${namestring}.root\ ${CMSSW_VERSION}\ ${USER}

cat > condor_${namestring}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_${namestring}.jdl
echo "Arguments = ${argument}" >> condor_${namestring}.jdl
cat >> condor_${namestring}.jdl << "EOF"
Output = simpleBH_$(Cluster)_$(Process).stdout
Error = simpleBH_$(Cluster)_$(Process).stderr
Log = simpleBH_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${namestring}.jdl
done
