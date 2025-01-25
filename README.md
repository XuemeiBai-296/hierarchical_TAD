###OnTAD
1.input: HiC-Pro "iced.matrix"; "abs.bed"
2.run 1-SparseToDense.sh
3.run 2-OnTAD.sh
4.output: ".tad"

###TH-score
1.input: ".tad"
2.run 3-TAD_domain.py calculate TAD domain
3.run 4-TH-score.py calculate TH-score

###ATAC-signal
1.input: ATAC-seq ".bw" files; TAD boundary file
2.run 5-ATAC_signal.sh
