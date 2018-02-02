#!/bin/bash
cd $P2X_ANALYSIS_CODE/Analysis
pwd
start=$1
end=$2
for i in $(seq $start $end)
do
file=/projects/prospect/data/AnalysisChallenge_PG4/Challenge_AD1_IBD/Run_${i}_DetSim.h5
if [ -e $file ] 
then
    echo "Processing Run_${i}_DetSim.h5"
    ./AnalyzeCalibrated  $file $DATA_CHALLENGE_OUTDIR/AD1_IBD/ AnalyzerConfig/IBDTreePlugin.cfg y
else
    echo "Run_${i}_DetSim.h5 ADI_IBD file does not exist"
fi


done