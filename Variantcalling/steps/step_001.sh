#!/bin/bash
input=$1
ref=$2
output=$3
mkdir -p $output
samples=$(ls $input | grep .bam$)
for i in $samples; do
    echo "___VARIANT_CALLING___"
    samtools index $input/$i
    gatk --java-options "-Xmx40g" HaplotypeCaller -R $ref -I $input/$i -O $output/${i}.vcf.gz
done



