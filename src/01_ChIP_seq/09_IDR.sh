for idx in MH_D_H3K9ac MH_N_H3K9ac MH_RE_H3K9ac
do
bsub -n 1 -q q2680v2 \
"idr \
-s ./01_narrowPeak/${idx}_1_peaks.narrowPeak ./01_narrowPeak/${idx}_2_peaks.narrowPeak \
--input-file-type narrowPeak \
--peak-list ./01_narrowPeak/${idx}_pooled_peaks.narrowPeak \
--rank signal.value \
--output-file ./04_IDR/${idx}.IDR.narrowPeak \
--output-file-type narrowPeak \
--plot \
-i 0.05 \
--peak-merge-method sum \
--use-best-multisummit-IDR"
sleep 2
done

for i in *.IDR.narrowPeak
do
idx=$(basename $i | sed 's/.IDR.narrowPeak//')
cut -f 1,2,3,4,5,6 $i | bedtools sort -i - | uniq > ../02_IDR_nonredundent_regions/${idx}.IDR.uniq.narrowPeak
done

bsub -n 1 -q normal \
"cat ./04_IDR/02_IDR_nonredundent_regions/*.narrowPeak |\
bedtools sort -i - | uniq | bedtools merge -d 0 -c 5 -o max > ./05_merged_IDR_peaks/MH_H3K9ac_unify.IDR.narrowPeak"



for i in DS NC
do
pool=${i}_bzip23_pooled_peaks.narrowPeak
rep1=${i}_bzip23_rep1_peaks.narrowPeak
rep2=${i}_bzip23_rep2_peaks.narrowPeak
bsub -n 1 -J idr -e idr.e -q q2680v2 \
"idr \
-s $rep1 $rep2 \
--peak-list $pool \
--input-file-type narrowPeak \
--output-file ${i}_bzip23_IDR.narrowPeak \
--output-file-type narrowPeak \
--plot \
--peak-merge-method avg \
--use-best-multisummit-IDR \
--idr-threshold 0.05"
sleep 2
done

for i in *_IDR.narrowPeak
do
idx=$(basename $i | sed 's/_IDR.narrowPeak//')
bsub -n 1 -q q2680v2 \
"cut -f 1,2,3 $i | bedtools sort -i - | uniq > ${idx}.uniq.bed"
sleep 2
done

cat DS_bzip23.uniq.bed NC_bzip23.uniq.bed | bedtools sort -i - | bedtools merge -i - | bedtools sort -i - > bzip23.unifiedPeaks.bed