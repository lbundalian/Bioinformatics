#!/bin/bash
input=$1
snpeff_util=$2
ref=$3
output=$4
mkdir -p $output

snps=$(ls $input | grep .vcf.gz$)

for snp in $snps; do
    echo "__COMMAND__"
    java -Xmx40g -jar $snpeff_util -v $ref $input/$snp > $output/${snp}.vcf
done
     



