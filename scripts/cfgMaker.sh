#!/usr/bin/sh

for i in 5 10 15 20 30 40 60 80 100
do
cat > simpleBH_E${i}Eta1p7.cfg << "EOF"
outFilePath = output/out_E${i}Eta1p7.root
EOF
echo "filePath = root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/"`ls /eos/uscms/store/user/\
chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/`"/0000" >> simpleBH_E${i}Eta1p7.cfg
cat >> simpleBH_E${i}Eta1p7.cfg << "EOF"
recoFileName = ntuples
nRuns = 80
deadfrac = 0.1
EOF
done
