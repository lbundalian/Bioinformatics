#!/bin/bash
input=$1
snpsift_util=$2
perl_util=$3
ref=$4
output=$5
mkdir -p $output

snps=$(ls $input | grep .vcf$)

for snp in $snps; do
    echo "__COMMAND__"
    cat $input/$snp | $perl_util | java -Xmx8g -jar $snpsift_util extractFields - "CHROM" "POS" "ID" "AF" "REF" "ALT" "ANN[*].GENE" "ANN[*].BIOTYPE" "ANN[*].IMPACT" "ANN[*].HGVS_P" "GEN[*].GT" "GEN[*].AD" "GEN[*].DP" "GEN[*].GQ" "GEN[*].PL" "DP" "QUAL"| awk '$8~/(HIGH|MODERATE)/ && $7=="protein_coding" && $9!="" && $15>9 && $16>29' - |  awk '!seen[$1,$2]++' - > $output/extracted_${snp}.txt
	
done
