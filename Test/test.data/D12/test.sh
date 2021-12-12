#!/bin/bash
read data
Dataset="$data"_mutated.fasta
Output="$data"_aligned.fasta
SPECIES_COUNT=$(grep ">" $Dataset | wc -l)


START=$(date +%s.%3N)
mafft  --auto --reorder --preservecase "$Dataset" > "$Output"
END=$(date +%s.%3N)
DIFF_ALIGNMENT=$(echo "$END - $START" | bc)


START=$(date +%s.%3N)
fasta=$Output
ORF=$(python FireMSA.py --input_file $fasta find_orf)
python FireMSA.py --input_file $fasta margin_trim --start $ORF --end None mask_sequence create_state_matrix 10 N remove_block_gaps 1 create_state_matrix 20 - remove_block_gaps 1 create_state_matrix 10 - create_state_matrix 5 N remove_block_gaps 3 find_frameshifts remove_seq_gaps 5 searc_stop_codon 3 save_msa FireMSA_Output terminate  ||  echo "Error encountered for custom script" 
END=$(date +%s.%3N)
DIFF_CUSTOM=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_CUSTOM=$(grep ">" FireMSA_Output.fasta | wc -l) || echo "Error encountered for custom script"


START=$(date +%s.%3N)
fasta=$Output
perl maxalign.pl -a "$fasta" || echo "Error encountered for maxalign script" 
END=$(date +%s.%3N)
DIFF_MAXALIGN=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_MAXAL=$(grep ">" heuristic.fsa | wc -l) || echo "Error encountered for maxalign script"


START=$(date +%s.%3N)
fasta=$Output
trimal -in $fasta -out Trimal_Output.fasta -automated1 || echo "Error encountered for trimal script"
END=$(date +%s.%3N)
DIFF_TRIMAL=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_TRIMAL=$(grep ">" Trimal_Output.fasta | wc -l) || echo "Error encountered for trimal script"


echo "Summary of Test" >> "results.txt"
echo "MAFFT ALignment: Time Elapsed $DIFF_ALIGNMENT sec" >> "results.txt"
echo "Custom Script : Time Elapsed $DIFF_CUSTOM sec Remaining Species $SP_CUSTOM" >> "results.txt"
echo "MaxALign Script : Time Elapsed $DIFF_MAXALIGN sec Remaining Species $SP_MAXAL" >> "results.txt"
echo "TrimAl Script : Time Elapsed $DIFF_TRIMAL sec Remaining Species $SP_TRIMAL" >> "results.txt"
