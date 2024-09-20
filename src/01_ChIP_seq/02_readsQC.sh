module load FastQC/0.11.9

for i in *_1.clean.fq.gz
do
idx=$(basename $i |sed 's/_1.clean.fq.gz//')
bsub -J QC_${idx} -n 1 -R span[hosts=1] -q q2680v2 \
"fastqc \
./${idx}_1.clean.fq.gz \
-o ./QC/

fastqc \
./${idx}_2.clean.fq.gz \
-o ./QC/"
sleep 2
done