bsub -n 1 -R span[hosts=1] -q q2680v2 \
"bowtie2-build MH63RS3.fa ./bowtieIdx/MH63"

for i in *_1.clean.fq.gz
do
idx=$(basename $i | sed 's/_1.clean.fq.gz//')
bsub -n 8 -e ../02_mapping/bt2log/${idx}.mapping_log -R span[hosts=1] -q q2680v2 \
"bowtie2 \
-p 8 \
-x  ~/genome_anno/mh63/bowtieIdx/MH63 \
-1 ${idx}_1.clean.fq.gz \
-2 ${idx}_2.clean.fq.gz \
--sensitive \
--no-unal \
-X 1000 \
| \
samtools view \
-@ 7 \
-f 2 \
-F 524 \
-b \
-o ../02_mapping/${idx}.bam"
sleep 2
done

for i in *.bam
do
idx=$(basename $i | sed 's/.bam//')
bsub -n 4 -R span[hosts=1] -q q2680v2 \
"samtools sort \
-l 9 \
-@ 3 \
-o ${idx}.sorted.bam \
${idx}.bam \
&& \
rm ${idx}.bam"
sleep 2
done