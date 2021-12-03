#!/bin/bash
START=$(date +%s)
perl translatorx_vLocal.pl -i d1_mutated.fasta
fasta="translatorx_res.nt_ali.fasta"
ORF=$(python FireMSA.py --input_file $fasta find_orf)
python FireMSA.py --input_file $fasta margin_trim --start $ORF --end None mask_sequence create_state_matrix 90 N remove_block_gaps 1 create_state_matrix 90 - remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 3 remove_seq_gaps 10 create_concensus match_concensus 20 save_msa output terminate
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "$DIFF sec"

