for i in early2cell late2cell 8cell ICM mESC
do
computeMatrix reference-point --referencePoint center -S ${i}.bw -R ./${i}/boundary_level1 ./${i}/boundary_level2 ./${i}/boundary_level3+ --beforeRegionStartLength 500000 --afterRegionStartLength 500000 --binSize 5000 --skipZeros -o ./${i}_matrix.mat.gz -p max

plotProfile -m ./${i}_matrix.mat.gz \
            -out ./${i}_Profile.pdf \
            #--plotTitle ${i}
done