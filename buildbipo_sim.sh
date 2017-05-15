#!/bin/bash
start=$1
end=$2
release=20161206
dataset=P50D

cd /projects/prospect/data/PG4_Simulation/P50D_Bi214


for i in $( ls *.root ) 
do
    root -b -q /home/jonesdc/prospect/macros/BuildBiPoTree.C+\(\"$i\",\"$num\",\"$dataset\",\"Phys_$release\"\)
done

