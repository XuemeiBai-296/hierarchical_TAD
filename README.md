### OnTAD
1. input: HiC-Pro "iced.matrix"; "abs.bed"
2. run 1-SparseToDense.sh
3. run 2-OnTAD.sh
4. output: ".tad"

### TH-score
1. input: ".tad"
2. run 3-TAD_domain.py calculate TAD domain
3. run 4-TH-score.py calculate TH-score

### ATAC-signal
1. input: ATAC-seq ".bw" files; TAD boundary file
2. run 5-ATAC_signal.sh

### the neural network
1.Calculation of boundary transition values for TAD in neighboring cell phases. eg. "late2cell-8cell_boundary_raise&reduce.txt"
2.Calculation of epigenetics signals overlapped with TAD boundary.
3."Xiewei_process_data.py" process the input data. output: "Xiewei_all_data.pkl", "Xiewei_all_data.pkl"
4. "model.py": First, encodding the high dimensional vectors for the current cell stage and the corresponding layer information into trainable hidden vectors of 16 dimensions using torch.nn.Embedding(). At the same time, concating the H3K4me3, H3K27me3, H3K9me3, and the next cell stage category, and construct 3 fully connected hidden layers through torch.nn.Linear(). The activation function of each layer is Tanh, and the outputs are the hidden vectors of 8 dimensions. Finally, the output implied vectors are predicted using a Sigmoid classifier, which outputs the probability of TAD boundary layer change for the next cell stage.
5. "xiewei_train.py"
