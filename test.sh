#!/bin/bash
start=$1
end=$2
release=20161206
dataset=P50D

for i in `seq $start $end`;
do 
    if [ $i -lt 10 ]; then
	num=00$i
    else
	num=0$i
    fi
    cd /home/prospect-collab/data/converted_data/Phys_$release/$dataset/series$num

    for i in $( ls *.root ) 
    do
	root -l -q /home/jonesdc/prospect/macros/test.C+\(\"$i\",\"$num\",\"$dataset\",\"Phys_$release\"\)
    done
done
