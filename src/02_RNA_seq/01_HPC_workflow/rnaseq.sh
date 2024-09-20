#!/bin/sh

step=$1

if [ -z $step ];then
echo "No step specified, input 'raw' or 'clean' to start with raw reads or clean reads respectively."
exit 1
fi

if [ $step != "raw" ] && [ $step != "clean" ];then
echo "Illegal step specified, input 'raw' or 'clean' to start with raw reads or clean reads respectively."
exit 1
fi

wd=$PWD
dir=$(cd $(dirname $0); pwd)

if [ $step == "raw" ];then
    if [ ! -d "00_raw" ]; then
    echo "Raw reads folder is not in the current directory, put raw reads files in '00_raw' and try again"
    exit 1
    else
    mkdir -p \
    $wd/01_clean \
    $wd/00_raw/trimlog \
    $wd/00_raw/QC

    #set reads quality value
    q=30

    cd $wd/00_raw
    for i in *_1.fq.gz
    do
    sh ${dir}/trim_qc.sh $i $q $wd &
    sleep 5
    done
    wait
    fi
elif [ $step == "clean" ];then
    if [ ! -d "01_clean" ]; then
    echo "Clean reads folder is not in the current directory, put clean reads files in '01_clean' and try again"
    exit 1
    fi
fi


#star_idx=/public/home/ychang/genome_anno/msu7/star_idx_150
#gff=/public/home/ychang/genome_anno/msu7/msu7.all.gff3
tar_idx=~/genome_anno/mh63/STAR_index
gff=~/genome_anno/mh63/MH63RS3.gff3

mkdir -p \
$wd/02_standard_DEG/01_mapping \
$wd/02_standard_DEG/02_stringtie \
$wd/02_standard_DEG/03_DEG

cd ${wd}/01_clean
for i in *_1.clean.fq.gz
do
sh ${dir}/star.sh ${i} ${wd} ${star_idx} ${gff} &
sleep 5
done

wait

cd ${wd}/02_standard_DEG/02_stringtie

sort count.txt -o count.txt
sort tabular.txt -o tabular.txt

bsub -K -J count -n 1 -R span[hosts=1]  -q q2680v2 \
"prepDE.py -i count.txt -g $wd/02_standard_DEG/03_DEG/gene_count.csv -l 150
rm transcript_count_matrix.csv"

bsub -K -J TPM -n 1 -R span[hosts=1]  -q q2680v2 \
"perl ${dir}/tab2TPM.pl -list tabular.txt -o $wd/02_standard_DEG/03_DEG/gene.tpm"

bsub -K -J TPM -n 1 -R span[hosts=1]  -q q2680v2 \
"perl ${dir}/tab2RPKM.pl -list tabular.txt -o $wd/02_standard_DEG/03_DEG/gene.fpkm"

cd $wd/02_standard_DEG/01_mapping

bsub -K -J bamsum -n 6 -R span[hosts=1]  -q q2680v2 \
"multiBamSummary bins \
-b *.sorted.bam \
-p max \
-bs 100 \
-o ./all_bam_summary.npz"

bsub -K -J bampca -n 1 -R span[hosts=1] -q q2680v2 \
"plotPCA \
-in all_bam_summary.npz \
--ntop 20000 \
-o ./bamsum_PCA.pdf"

bsub -K -J bamcorr -n 1 -R span[hosts=1] -q q2680v2 \
"plotCorrelation \
-in all_bam_summary.npz \
--corMethod pearson \
--skipZeros \
--plotTitle 'Pearson Correlation of Average Scores Per Transcript' \
--whatToPlot scatterplot \
--removeOutliers \
-o bamcorr_scatterplot.pdf"

bsub -K -J bamcorr -n 1 -R span[hosts=1] -q q2680v2 \
"plotCorrelation \
-in all_bam_summary.npz \
--corMethod pearson \
--skipZeros \
--plotTitle 'Pearson Correlation of Average Scores Per Transcript' \
--whatToPlot heatmap \
--plotNumbers \
--removeOutliers \
-o bamcorr_heatmap.pdf \
--outFileCorMatrix bamcorr_Pcc.tab"

ln -s $wd/sample_info/*.csv ${wd}/02_standard_DEG/03_DEG

cd ${wd}/02_standard_DEG/03_DEG

#bsub -K -J Deseq2 -n 4 -R span[hosts=1] -q normal \
#"source activate r
#Rscript ${dir}/deseq2.R"
