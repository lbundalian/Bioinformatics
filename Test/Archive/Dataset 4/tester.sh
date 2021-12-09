#!/bin/bash
Dataset=d4_mutated.fasta
SPECIES_COUNT=$(grep ">" $Dataset | wc -l)

START=$(date +%s)
perl translatorx_vLocal.pl -i $Dataset
fasta="translatorx_res.nt_ali.fasta"
ORF=$(python FireMSA.py --input_file $fasta find_orf)
python FireMSA.py --input_file $fasta margin_trim --start $ORF --end None mask_sequence create_state_matrix 10 N remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 1 create_state_matrix 10 - remove_block_gaps 3 remove_seq_gaps 10 find_frameshifts create_concensus match_concensus 20 search_stop_codon 3 save_msa output terminate
END=$(date +%s)
DIFF_CUSTOM=$(echo "$END - $START" | bc)
SP_CUSTOM=$(grep ">" output.fasta | wc -l)
#echo "Custom Script : $DIFF sec"


#!/bin/bash
START=$(date +%s)
perl translatorx_vLocal.pl -i $Dataset
fasta="translatorx_res.nt_ali.fasta"
perl maxalign.pl -a $fasta
END=$(date +%s)
DIFF_MAXALIGN=$(echo "$END - $START" | bc)
SP_MAXAL=$(grep ">" heuristic.fsa | wc -l)
#echo "MaxAlign Script : $DIFF sec"


#!/bin/bash
START=$(date +%s)
perl translatorx_vLocal.pl -i $Dataset
fasta="translatorx_res.nt_ali.fasta"
trimal -in translatorx_res.nt_ali.fasta -out trimal_output.fasta -automated1
END=$(date +%s)
DIFF_TRIMAL=$(echo "$END - $START" | bc)
SP_TRIMAL=$(grep ">" trimal_output.fasta | wc -l)
#echo "MaxAlign Script : $DIFF sec"


echo "Summary of Test"
echo "Custom Script : Time Elapsed $DIFF_CUSTOM sec Remaining Species $SP_CUSTOM"
echo "MaxALign Script : Time Elapsed $DIFF_MAXALIGN sec Remaining Species $SP_MAXAL"
echo "TrimAl Script : Time Elapsed $DIFF_TRIMAL sec Remaining Species $SP_TRIMAL"
