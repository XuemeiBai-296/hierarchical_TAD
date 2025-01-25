###TH score
import os

file_path = "./GSE82185_17_Nature_Xiewei/20kb/"
dirs = os.listdir(file_path)

for s in dirs:
    os.system("bedtools coverage -mean -a gene_body.bed -b ./GSE82185_17_Nature_Xiewei/20kb/"+s+"/domain.bed > ./"+s+"/TH_score")