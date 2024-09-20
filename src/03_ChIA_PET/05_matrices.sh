ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*pool*/*ipet ./
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*pool*/*spet ./

ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*rep1*/*ipet ./
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*rep1*/*spet ./

ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*rep2*/*ipet ./
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/*rep2*/*spet ./

for i in D N RE
do
for j in rep1 rep2 pooled
do
bsub -n 1 -q normal \
"
cat MH_H3K9ac_${i}_${j}.ipet MH_H3K9ac_${i}_${j}.spet > MH_H3K9ac_${i}_${j}.sipet
"
sleep 2
done
done

for i in MH_H3K9ac_D MH_H3K9ac_N MH_H3K9ac_RE
do
for bs in 1 5 10 50 100
do
bsub -n 1 -q q2680v2 -e bedpe2Matrix.e -R span[hosts=1] \
"bedpe2Matrix \
--binsize ${bs}000 \
--chrsizes ~/genome_anno/mh63/MH63RS3.chrlength.txt \
--ifile $i.sipet \
--oprefix ../02_Hic_files/$i.${bs}k \
--all \
--binoffset 1 \
--matrix-format complete"
sleep 2
done
done

conda activate chiapet

for i in  MH_H3K9ac_D MH_H3K9ac_N MH_H3K9ac_RE
do
bsub -n 1 -e e.e  -q q2680v2 \
"hicConvertFormat \
-m ${i}*matrix \
-o ../03_H5_cool/${i}.mcool \
--inputFormat hicpro \
--outputFormat mcool \
-bf ${i}*abs.bed"
sleep 2
done

prefix1=MH_H3K9ac_N
prefix2=MH_H3K9ac_D
prefix3=MH_H3K9ac_RE
prefix4=bzip23_pooled
outdir=../04_normalized_matrix
for i in 1 5 10 50 100
do
bsub -n 1 -e e.e -q q2680v2 \
"hicNormalize \
-m ${prefix1}.mcool::/resolutions/${i}000 \
${prefix2}.mcool::/resolutions/${i}000 \
${prefix3}.mcool::/resolutions/${i}000 \
${prefix4}.mcool::/resolutions/${i}000 \
-n smallest \
-o $outdir/${prefix1}.${i}kb.cool \
$outdir/${prefix2}.${i}kb.cool \
$outdir/${prefix3}.${i}kb.cool \
$outdir/${prefix4}.${i}kb.cool"
sleep 2
done

for i in  MH_H3K9ac_D MH_H3K9ac_N MH_H3K9ac_RE bzip23_pooled
do
bsub -n 1 -e e.e  -q q2680v2 \
"hicConvertFormat \
-m ${i}*.cool \
-o ./${i}.mcool \
--inputFormat cool \
--outputFormat mcool"
sleep 2
done

for bs in 1 5 10 50 100
do
bsub -n 1 -q q2680v2 -e mtxcomp.e \
"hicCompareMatrices -m \
./bzip23_pooled.mcool::/resolutions/${bs}000 \
./MH_H3K9ac_D.mcool::/resolutions/${bs}000 \
--operation log2ratio \
-o ../05_matrixCompare/bzip23_vs_MH_DS_H3K9ac.matrixCompare.${bs}000.log2ratio.cool"
sleep 2
done

for bs in 1 5 10 50 100
do
bsub -n 1 -q q2680v2 -e mtxcomp.e \
"hicCompareMatrices -m \
./MH_H3K9ac_D.mcool::/resolutions/${bs}000 \
./MH_H3K9ac_N.mcool::/resolutions/${bs}000 \
--operation log2ratio \
-o ../05_matrixCompare/MH_H3K9ac_DS_vs_NC.matrixCompare.${bs}000.log2ratio.cool"
sleep 2
done

for bs in 1 5 10 50 100
do
bsub -n 1 -q q2680v2 -e mtxcomp.e \
"hicCompareMatrices -m \
./MH_H3K9ac_RE.mcool::/resolutions/${bs}000 \
./MH_H3K9ac_D.mcool::/resolutions/${bs}000 \
--operation log2ratio \
-o ../05_matrixCompare/MH_H3K9ac_RE_vs_DS.matrixCompare.${bs}000.log2ratio.cool"
sleep 2
done

for bs in 1 5 10 50 100
do
bsub -n 1 -q q2680v2 -e mtxcomp.e \
"hicCompareMatrices -m \
./MH_H3K9ac_RE.mcool::/resolutions/${bs}000 \
./MH_H3K9ac_N.mcool::/resolutions/${bs}000 \
--operation log2ratio \
-o ../05_matrixCompare/MH_H3K9ac_RE_vs_NC.matrixCompare.${bs}000.log2ratio.cool"
sleep 2
done

for i in  MH_H3K9ac_DS_vs_NC MH_H3K9ac_RE_vs_DS MH_H3K9ac_RE_vs_NC bzip23_vs_MH_DS_H3K9ac
do
bsub -n 1 -e e.e  -q q2680v2 \
"hicConvertFormat \
-m ${i}*.cool \
-o ./${i}.matrixCompare.log2ratio.mcool \
--inputFormat cool \
--outputFormat mcool"
sleep 2
done

bsub -n 1 -q q2680v2 -e e.e \
"hicPlotMatrix \
--vMax 3 \
--vMin -3 \
--colorMap RdBu_r \
-m ./bzip23_vs_MH_DS_H3K9ac.matrixCompare.log2ratio.mcool::/resolutions/100000 \
-o ../06_matrix_view/GW_bzip23_vs_MH_DS_H3K9ac.matrixCompare.100000.log2ratio.heatmap.pdf \
--dpi 600"

bsub -n 1 -q q2680v2 -e e.e \
"hicPlotMatrix \
--vMax 3 \
--vMin -3 \
--colorMap RdBu_r \
-m ./MH_H3K9ac_DS_vs_NC.matrixCompare.log2ratio.mcool::/resolutions/100000 \
-o ../06_matrix_view/GW_MH_H3K9ac_DS_vs_NC.matrixCompare.100000.log2ratio.heatmap.pdf \
--dpi 600"

bsub -n 1 -q q2680v2 -e e.e \
"hicPlotMatrix \
--vMax 3 \
--vMin -3 \
--colorMap RdBu_r \
-m ./MH_H3K9ac_RE_vs_NC.matrixCompare.log2ratio.mcool::/resolutions/100000  \
-o ../06_matrix_view/GW_MH_H3K9ac_RE_vs_NC.matrixCompare.100000.log2ratio.heatmap.pdf \
--dpi 600"

bsub -n 1 -q q2680v2 -e e.e \
"hicPlotMatrix \
--vMax 3 \
--vMin -3 \
--colorMap RdBu_r \
-m ./MH_H3K9ac_RE_vs_DS.matrixCompare.log2ratio.mcool::/resolutions/100000  \
-o ../06_matrix_view/GW_MH_H3K9ac_RE_vs_DS.matrixCompare.100000.log2ratio.heatmap.pdf \
--dpi 600"

for i in D N RE
do
for j in rep1 rep2
do
bsub -n 1 -q normal \
"
cat MH_H3K9ac_${i}_${j}.ipet MH_H3K9ac_${i}_${j}.spet > MH_H3K9ac_${i}_${j}.sipet
"
sleep 2
done
done

for i in MH_H3K9ac_D MH_H3K9ac_N MH_H3K9ac_RE
do
for j in rep1 rep2
do
for bs in 5 50 100
do
bsub -n 1 -q q2680v2 -e bedpe2Matrix.e -R span[hosts=1] \
"bedpe2Matrix \
--binsize ${bs}000 \
--chrsizes ~/genome_anno/mh63/MH63RS3.chrlength.txt \
--ifile ${i}_${j}.sipet \
--oprefix ../02_mcool/${i}_${j}.${bs}k \
--all \
--binoffset 1 \
--matrix-format complete"
sleep 2
done
done
done

conda activate chiapet

for i in MH_H3K9ac_D_rep1 MH_H3K9ac_N_rep1 MH_H3K9ac_RE_rep1 MH_H3K9ac_D_rep2 MH_H3K9ac_N_rep2 MH_H3K9ac_RE_rep2
do
bsub -n 1 -e e.e  -q q2680v2 \
"hicConvertFormat \
-m ${i}*matrix \
-o ./${i}.mcool \
--inputFormat hicpro \
--outputFormat mcool \
-bf ${i}*abs.bed"
sleep 2
done

prefix1=MH_H3K9ac_N
prefix2=MH_H3K9ac_D
prefix3=MH_H3K9ac_RE
outdir=../03_normalize
for i in 5 50 100
do
bsub -n 1 -e e.e -q q2680v2 \
"hicNormalize \
-m ${prefix1}_rep1.mcool::/resolutions/${i}000 \
${prefix1}_rep2.mcool::/resolutions/${i}000 \
${prefix2}_rep1.mcool::/resolutions/${i}000 \
${prefix2}_rep2.mcool::/resolutions/${i}000 \
${prefix3}_rep1.mcool::/resolutions/${i}000 \
${prefix3}_rep2.mcool::/resolutions/${i}000 \
-n smallest \
-o $outdir/${prefix1}_rep1.${i}kb.cool \
$outdir/${prefix1}_rep2.${i}kb.cool \
$outdir/${prefix2}_rep1.${i}kb.cool \
$outdir/${prefix2}_rep2.${i}kb.cool \
$outdir/${prefix3}_rep1.${i}kb.cool \
$outdir/${prefix3}_rep2.${i}kb.cool"
sleep 2
done

for i in  MH_H3K9ac_D MH_H3K9ac_N MH_H3K9ac_RE
do
for j in rep1 rep2
do
bsub -n 1 -e e.e  -q q2680v2 \
"hicConvertFormat \
-m ${i}_${j}*.cool \
-o ./${i}_${j}.norm.mcool \
--inputFormat cool \
--outputFormat mcool"
sleep 2
done
done

for i in D N RE
do
for j in 5 50 100
do
bsub -n 1 -q q2680v2 \
"
hicCorrelate \
-m MH_H3K9ac_${i}_rep1.norm.mcool::/resolutions/${j}000 MH_H3K9ac_${i}_rep2.norm.mcool::/resolutions/${j}000 \
--log1p \
--plotNumbers \
--range 8000:2000000 \
--method pearson \
-oh ../04_correlation/MH_H3K9ac_${i}.res${j}k.PCC.heatmap.pdf \
-os ../04_correlation/MH_H3K9ac_${i}.res${j}k.PCC.scatter.pdf
"
sleep 2
done
done
