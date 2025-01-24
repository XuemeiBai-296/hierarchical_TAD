#!bin/bash
r="20000"
Res="20kb"
dirs="GSE82185_17_Nature_Xiewei"
out_dir="./1-SparseToDense/dence_matrix/${dirs}/${Res}"
cd $out_dir
python ~/HiC-Pro_3.1.0/bin/utils/sparseToDense.py -b  ./early2cell_20000_abs.bed ./early2cell_20000_iced.matrix -o ${i}_${Res}_dence.matrix --perchr
cd ..
done

