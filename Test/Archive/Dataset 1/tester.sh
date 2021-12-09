#!/bin/bash
Dataset=d1_mutated.fasta
Output=d1_aligned.fasta
SPECIES_COUNT=$(grep ">" $Dataset | wc -l)


START=$(date +%s.%3N)
mafft  --auto --reorder --preservecase "$Dataset" > "$Output"
END=$(date +%s.%3N)
DIFF_ALIGNMENT=$(echo "$END - $START" | bc)


START=$(date +%s.%3N)
fasta=$Output
ORF=$(python FireMSA.py --input_file $fasta find_orf)
python FireMSA.py --input_file $fasta margin_trim --start $ORF --end None mask_sequence create_state_matrix 10 N remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 3 remove_seq_gaps 10 find_frameshifts create_concensus match_concensus 20 search_stop_codon 3 save_msa FireMSA_Output terminate
END=$(date +%s.%3N)
DIFF_CUSTOM=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_CUSTOM=$(grep ">" FireMSA_Output.fasta | wc -l)


START=$(date +%s.%3N)
fasta=$Output
perl maxalign.pl -a "$fasta"
END=$(date +%s.%3N)
DIFF_MAXALIGN=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_MAXAL=$(grep ">" heuristic.fsa | wc -l)


START=$(date +%s.%3N)
fasta=$Output
trimal -in $fasta -out Trimal_Output.fasta -automated1
END=$(date +%s.%3N)
DIFF_TRIMAL=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_TRIMAL=$(grep ">" Trimal_Output.fasta | wc -l)


echo "Summary of Test"
echo "MAFFT ALignment: Time Elapsed $DIFF_ALIGNMENT sec"
echo "Custom Script : Time Elapsed $DIFF_CUSTOM sec Remaining Species $SP_CUSTOM"
echo "MaxALign Script : Time Elapsed $DIFF_MAXALIGN sec Remaining Species $SP_MAXAL"
echo "TrimAl Script : Time Elapsed $DIFF_TRIMAL sec Remaining Species $SP_TRIMAL"
