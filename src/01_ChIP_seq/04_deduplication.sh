for i in *.sorted.bam
do
idx=$(basename $i | sed 's/.sorted.bam//')
bsub -n 1 -R span[hosts=1] -q q2680v2 \
"picard MarkDuplicates \
I=${idx}.sorted.bam \
O=${idx}.rmdup.bam \
M=${idx}.dupmetrics.txt \
REMOVE_DUPLICATES=true && \
rm ${idx}.sorted.bam"
sleep 2
done