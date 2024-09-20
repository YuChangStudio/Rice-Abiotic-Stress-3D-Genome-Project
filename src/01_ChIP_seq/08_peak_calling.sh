for i in *H3K9ac*IP.rmdup.bam
do
idx=$(basename $i | sed 's/_IP.rmdup.bam//')
bsub  -n 1 -R span[hosts=1] -q q2680v2 \
"macs2 callpeak \
--call-summits \
-t ${idx}_IP.rmdup.bam \
-c MH_input.rmdup.bam \
-f BAMPE \
-g 3.9e8 \
-p 0.01 \
--keep-dup all \
--max-gap 500 \
-n ${idx} \
--outdir ../04_peakcall_summit"
sleep 2
done

for i in *H3K9ac*IP.sorted.pooled.bam
do
idx=$(basename $i | sed 's/_IP.sorted.pooled.bam//')
bsub  -n 1 -R span[hosts=1] -q q2680v2 \
"macs2 callpeak \
--call-summits \
-t ${idx}_IP.sorted.pooled.bam \
-c MH_input.rmdup.bam \
-f BAMPE \
-g 3.9e8 \
-p 0.01 \
--keep-dup all \
-n ${idx}_pooled \
--max-gap 500 \
--outdir ../04_peakcall_summit"
sleep 2
done


for i in DS NC
do
for j in pooled rep1 rep2
do
ip=${i}_bzip23IP_${j}.rmdup.bam
input=${i}_input.rmdup.bam
out=${i}_bzip23_${j}
bsub  -n 1 -R span[hosts=1] -e e.e -q q2680v2 \
"macs2 callpeak \
--call-summits \
-t $ip \
-c $input \
-f BAMPE \
-g 3.9e8 \
-p 0.05 \
--keep-dup all \
-n $out \
--outdir ../05_peakcall"
sleep 2
done
done
