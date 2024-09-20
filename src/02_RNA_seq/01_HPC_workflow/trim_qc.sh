#!/bin/sh
i=$1
q=$2
wd=$3
idx=$(basename $i |sed 's/_1.fq.gz//')

bsub -K -J trim_${idx} -n 4 -e ./trimlog/${idx}.trimlog -R span[hosts=1] -q normal \
"fastp \
-q ${q} \
-i ${idx}_1.fq.gz \
-o ${wd}/01_clean/${idx}_1.clean.fq.gz \
-I ${idx}_2.fq.gz \
-O ${wd}/01_clean/${idx}_2.clean.fq.gz \
-w 4 \
-j ./trimlog/${idx}.josn \
-h ./trimlog/${idx}.html \
-R ${idx}"


bsub -K -J QC_${idx} -n 1 -R span[hosts=1] -q normal \
"fastqc \
${wd}/01_clean/${idx}_1.clean.fq.gz \
-o ./QC/

fastqc \
${wd}/01_clean/${idx}_2.clean.fq.gz \
-o ./QC/"
