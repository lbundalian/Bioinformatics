#!/bin/bash
for i in {3..12}
do
 cp AlignmentScore.py D$i
 cd D$i
 echo "Alignment Scores" >> Scores.txt
 echo "Original aligned file" >> Scores.txt
 python AlignmentScore.py -i d${i}_aligned.fasta >> Scores.txt
 echo "Custom" >> Scores.txt
 python AlignmentScore.py -i FireMSA_Output.fasta >> Scores.txt
 echo "Trimal" >> Scores.txt
 python AlignmentScore.py -i Trimal_Output.fasta >> Scores.txt
 echo "Maxal" >> Scores.txt
 python AlignmentScore.py -i heuristic.fsa >> Scores.txt
 echo "End"
 cd ..
done
