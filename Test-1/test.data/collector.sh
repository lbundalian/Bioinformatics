#!/bin/bash
for i in {1..12}
do
cd Datasets/D$i/
cp d${i}_original.fasta ~/All/
cd ..
cd ..
done
