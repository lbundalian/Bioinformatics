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
fasta=$Dataset
ORF=$(python FireMSA.py --input_file $fasta find_orf)
python FireMSA.py --input_file $fasta margin_trim --start $ORF --end None search_stop_codon 3 save_msa FireMSA_Output terminate ||  echo "Error encountered for custom script" 
END=$(date +%s.%3N)
DIFF_CUSTOM=$(echo "$END - $START + $DIFF_ALIGNMENT" | bc)
SP_CUSTOM=$(grep ">" FireMSA_Output.fasta | wc -l) || echo "Error encountered for custom script"

echo "Summary of Test" >> "results.txt"
echo "MAFFT ALignment: Time Elapsed $DIFF_ALIGNMENT sec" >> "results.txt"
echo "Custom Script : Time Elapsed $DIFF_CUSTOM sec Remaining Species $SP_CUSTOM" >> "results.txt"

