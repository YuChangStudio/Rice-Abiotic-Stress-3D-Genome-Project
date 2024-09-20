bsub -n 4 -R span[hosts=1] -q normal \
"computeMatrix \
scale-regions \
-p 4 \
--scoreFileName *.bw  \
--regionsFileName *.bed \
--startLabel TSS \
--endLabel TES \
--upstream 2000 \
--downstream 2000 \
--outFileName H3K9ac_on_all.SR.matrix.gz \
--regionBodyLength 2000 \
--binSize 10"

bsub -n 1 -R span[hosts=1] -q normal \
"plotHeatmap \
--matrixFile H3K9ac_on_all.SR.matrix.gz  \
--outFileName H3K9ac_on_all.SR.pdf \
--alpha 0.7 \
--plotTitle H3K9ac_on_all_genes \
--colorMap plasma_r"

bsub -n 1 -R span[hosts=1] -q normal \
"plotHeatmap \
--matrixFile H3K9ac_on_all.SR.matrix.gz  \
--outFileName H3K9ac_on_all.SR.pdf \
--alpha 0.7 \
--plotTitle H3K9ac_on_all_genes \
--colorMap plasma_r"



for i in {1..7}
do
incluster=Fig3.SET${i}.bedpe
out=Fig3.SET${i}.loopspan.bed
awk -F "\t" '{if($1==$4){OFS="\t";print $1,int(($3+$2)/2),int(($6+$5)/2)}}' $incluster > $out
done

for i in {1..3}
do
incluster=Fig5.SET${i}.bedpe
out=Fig5.SET${i}.loopspan.bed
awk -F "\t" '{if($1==$4){OFS="\t";print $1,int(($3+$2)/2),int(($6+$5)/2)}}' $incluster > $out
done

for i in DSs nDSs
do
incluster=Fig5.SET1_${i}.bedpe
out=Fig5.SET1_${i}.loopspan.bed
awk -F "\t" '{if($1==$4){OFS="\t";print $1,int(($3+$2)/2),int(($6+$5)/2)}}' $incluster > $out
done

bsub -n 4 -R span[hosts=1] -q q2680v2 \
"computeMatrix \
scale-regions \
-p 4 \
--scoreFileName *.bw  \
--regionsFileName *.bed \
--startLabel TSS \
--endLabel TES \
--upstream 1000 \
--downstream 1000 \
--outFileName bzip23.onSETs.SR10bp.matrix.gz \
--regionBodyLength 10000 \
--binSize 10 && \
plotHeatmap \
--matrixFile bzip23.onSETs.SR10bp.matrix.gz  \
--outFileName bzip23.onSETs.SR10bp.heatmap.pdf \
--plotType se \
--alpha 1 \
--perGroup \
--colorMap inferno_r \
--plotTitle bzip23.onSETs.SR10bp"