#!/bin/sh
i=$1
wd=$2
star_idx=$3
gff=$4

idx=$(basename $i |sed 's/_1.clean.fq.gz//')

bsub -K -J star_${idx} -n 8 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/01_clean
STAR \
--runThreadN 8 \
--genomeDir ${star_idx} \
--readFilesIn ${idx}_1.clean.fq.gz ${idx}_2.clean.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${wd}/02_standard_DEG/01_mapping/${idx}."

bsub -K -J sam2bam_${idx} -n 2 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/02_standard_DEG/01_mapping
samtools view \
-b \
-f 2 \
-F 524 \
-q 5 \
-@ 1 \
-o ${idx}.bam \
${idx}.Aligned.out.sam"

bsub -K -J sortbam_${idx} -n 2 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/02_standard_DEG/01_mapping
samtools sort \
-l 9 \
-@ 1 \
-o ${idx}.sorted.bam \
${idx}.bam
rm ${idx}.Aligned.out.sam ${idx}.bam"

bsub  -K -J bamIDX_${idx} -n 1 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/02_standard_DEG/01_mapping
samtools index \
-b ${idx}.sorted.bam \
${idx}.sorted.bam.bai"

bsub -K -J bw_${idx} -n 8 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/02_standard_DEG/01_mapping
bamCoverage \
--numberOfProcessors 8 \
--normalizeUsing RPKM \
-b ${idx}.sorted.bam \
-o ${idx}.RPKM.bw \
-bs 1 \
--smoothLength 3"

bsub -K -J stringtie_${idx} -n 2 -R span[hosts=1] -q q2680v2 \
"cd ${wd}/02_standard_DEG/02_stringtie
stringtie \
--rf \
-p 2 \
-c 2 \
-e \
-G ${gff} \
-o ${idx}.gtf \
-A ${idx}.tab  \
${wd}/02_standard_DEG/01_mapping/${idx}.sorted.bam

echo ${idx} ${wd}/02_standard_DEG/02_stringtie/${idx}.gtf >> count.txt
echo ${idx} ${wd}/02_standard_DEG/02_stringtie/${idx}.tab >> tabular.txt"
