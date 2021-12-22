#!/bin/bash
for i in {1..12}
do
 echo "D$i" >> Consensus.txt
 python FireConsensus.py --seed d${i}_guidance.fasta --fasta d${i}_original.fasta score_consensus >> Consensus.txt
done
