#!/usr/bin/sh

for df in 01 03 05 07
do
for i in 5 10 15 20 30 40 60 80 100
do
echo "outFilePath = out_E${i}Eta1p7_df${df}.root" > simpleBH_E${i}Eta1p7_df${df}.cfg
echo "filePath = root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/"`ls /eos/uscms/store/user/\
chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/`"/0000" >> simpleBH_E${i}Eta1p7_df${df}.cfg
cat >> simpleBH_E${i}Eta1p7_df${df}.cfg << "EOF"
recoFileName = ntuples
nRuns = 80
EOF
echo "deadfrac = 0.${df}" >> simpleBH_E${i}Eta1p7_df${df}.cfg
done
done
