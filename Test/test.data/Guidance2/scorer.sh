#!/bin/bash
for i in {1..12}
do
echo D$i >> Scores.txt
python AlignmentScore.py -i d${i}_guidance.fasta >> Scores.txt
echo "Done D$i"
done
