#!/bin/bash
###1.define input
dirs="GSE82185_17_Nature_Xiewei"
res="20kb"
r="20000"
min="2"
max="100"
chr_int=./chrom_mm10.txt
matrix_dr=./1-SparseToDense/dence_matrix/${dirs}/${res}/
out_dr=./2-OnTAD_result/${dirs}/${res}/
cd ${out_dr}
for ii in $( ls ${matrix_dr})
do
/mnt/f/Early_embryo_TAD_hierarchy/Code/2_OnTAD/OnTAD-master/OnTAD ${matrix_dr}${ii} -penalty 0.1 -minsz ${min} -maxsz ${max} -o ${ii}_pen0.1_max${max}_${chr} -bedout ${chr} ${chr_size} ${r}
done 
