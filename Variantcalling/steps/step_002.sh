#!/bin/bash
input=$1
ref=$2
var=$3
output1=$4
mkdir -p $output1
vcfs=$(ls $input | grep .vcf.gz$)

for vc in $vcfs; do
    echo "__GET_VARIANTS__"
    gatk --java-options "-Xmx40g" SelectVariants -R $ref -V $input/$vc --select-type-to-include $var -O $output1/${vc}
done
     



