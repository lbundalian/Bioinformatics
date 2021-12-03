#!/bin/bash
START=$(date +%s)
perl translatorx_vLocal.pl -i d1_mutated.fasta
fasta="translatorx_res.nt_ali.fasta"
ORF=$(python FireMSA.py --input_file $fasta find_orf)
echo $fasta
$(python FireMSA.py --input_file $fasta.fasta margin_trim --start $ORF --end None mask_sequence create_state_matrix 90 N remove_block_gaps 1 create_state_matrix 90 - remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 3 remove_seq_gaps 25 save_msa output)



END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."

