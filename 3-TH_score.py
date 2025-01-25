###TAD domain
import numpy as np
import os

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')
    
chrname = []
for i in range(1,20):
    chrname.append(str(i))
chrname.append("X")

file_path = "./GSE82185_17_Nature_Xiewei/20kb/"
dirs = os.listdir(file_path)
# for s in dirs:
    # os.system("mkdir -p 5-OnTAD_domain/"+s)

for s in dirs:
    domain = []
    for chrom in chrname:
        f = open("./GSE82185_17_Nature_Xiewei/20kb/"+s+"/chr"+chrom+"_"+s+"_20kb_dence.matrix_pen0.1_max100_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            domain.append(["chr"+chrom,str((int(sep_line[0])-1)*20000),str((int(sep_line[1])-1)*20000)])
        f.close()
        
    savetxt("./GSE82185_17_Nature_Xiewei/20kb/"+s+"/domain.bed",domain)

###TH score
import os

file_path = "./GSE82185_17_Nature_Xiewei/20kb/"
dirs = os.listdir(file_path)

# for s in dirs:
#     os.system("mkdir -p"+s)

for s in dirs:
    os.system("bedtools coverage -mean -a gene_body.bed -b ./GSE82185_17_Nature_Xiewei/20kb/"+s+"/domain.bed > ./"+s+"/TH_score")