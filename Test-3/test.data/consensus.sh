#!/bin/bash
for i in {2..5}
do
 cp FireConsensus.py D$i
 cd D$i
 echo "Consensus Scores" >> Consensus.txt
 echo "Custom" >> Consensus.txt
 python FireConsensus.py --seed FireMSA_Output.fasta --fasta d${i}_original.fasta score_consensus >> Consensus.txt
 echo "Trimal" >> Consensus.txt
 python FireConsensus.py --seed Trimal_Output.fasta --fasta d${i}_original.fasta score_consensus >> Consensus.txt
 echo "Maxal" >> Consensus.txt
 python FireConsensus.py -seed heuristic.fsa --fasta d${i}_original.fasta score_consensus >> Consensus.txt
 echo "D$i done"
 cd ..
done
