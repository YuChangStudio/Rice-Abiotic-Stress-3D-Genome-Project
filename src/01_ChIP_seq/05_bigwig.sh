for i in *.sorted.pooled.bam
do
idx=$(basename $i |sed 's/.sorted.pooled.bam//')
bsub -n 8 -R span[hosts=1] -q q2680v2 \
"bamCoverage \
--numberOfProcessors 8 \
--normalizeUsing RPGC \
-b $i \
-o $idx.bw \
-bs 1 \
--effectiveGenomeSize 395765488"
sleep 2
done