for i in *_1.fq.gz
do
idx=$(basename $i | sed 's/_1.fq.gz//')
bsub -J trim_${idx} -n 4  -e ./trimlog/${idx}.trimlog -R span[hosts=1]  -q q2680v2 \
"fastp \
-q 15 \
-i ${idx}_1.fq.gz \
-o ../01_clean/${idx}_1.clean.fq.gz \
-I ${idx}_2.fq.gz \
-O ../01_clean/${idx}_2.clean.fq.gz \
-w 4 \
-j ./trimlog/${idx}.josn \
-h ./trimlog/${idx}.html \
-R ${idx}"
sleep 2
done