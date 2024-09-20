for i in *
do
cd $i
bsub -n 6 -R span[hosts=1] -q q2680v2 \
"samtools merge -@ 5 -l 9 - *.sam | \
samtools view -@ 5 -b - | \
samtools sort -@ 5 -l 9 -o $i.sorted.bam && \
samtools index -@ 5 $i.sorted.bam $i.sorted.bam.bai && \
rm *.sam
"
cd ..
sleep 2
done


for i in *
do
cd $i
bsub -n 8 -R span[hosts=1] -q q2680v2 \
"bamCoverage \
--numberOfProcessors 8 \
--normalizeUsing RPGC \
-b $i.sorted.bam \
-o $i.bw \
-bs 10 \
--effectiveGenomeSize 395765488 \
--minFragmentLength 18 \
--maxFragmentLength 1000"
sleep 2
cd ..
done