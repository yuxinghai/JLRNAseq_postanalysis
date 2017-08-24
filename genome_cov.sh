#!/usr/bin/bash
#KD_SS_PSF-v-Control.tab is primary result of IRfinder
NAMES=(psf_IR_L10k psf_IR_L20k psf_IR_L30k psf_IR_L40k psf_IR_L50k psf_IR_L60k psf_IR_new  KD_SS_PSF-v-Control.tab)
numSamples=${#NAMES[@]}
for (( i=0; i<${numSamples}; i++ )); do
    name=${NAMES[$i]}
    python gene_position.py -i ${name} -r hg19_37_gene.gtf -o ${name}_genomeCov.bed
    sed -r -i 's/chr[0-9]+/chrN/' ${name}_genomeCov.bed
    awk -v OFS="\t" '{print $1,$2,$3}' ${name}_genomeCov.bed >${name}_genomeCov_I1.bed 
    bedtools coverage -a my_genome.bed -b ${name}_genomeCov_I1.bed -d >${name}_pre_mygenome.txt
done

