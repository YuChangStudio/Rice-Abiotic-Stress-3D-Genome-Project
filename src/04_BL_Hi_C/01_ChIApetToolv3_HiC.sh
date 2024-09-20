conda activate chiapet

REF=~/genome_anno/mh63/bwaid/MH63RS3bwaid
chrlength=~/genome_anno/mh63/MH63RS3.chrlength.txt
cytoband=~/genome_anno/mh63/MH63RS3.cytoBandIdeo.txt
gs=3.9E8
OUT=../03_ChIApetTool_out

for i in *_1.fq
do
idx=$(basename $i |sed 's/_1.fq//')
bsub -n 8 -q smp -o cptv3.o -e cptv3.e -R "rusage[mem=64GB]" -R span[hosts=1] \
"java -jar ~/software/ChIA-PET_Tool_V3-master/ChIA-PET.jar \
--mode 1 \
--thread 8 \
--start_step 1 \
--fastq1 ${idx}_1.fq \
--fastq2 ${idx}_2.fq \
--minimum_linker_alignment_score 19 \
--linker ../BL_MboI.txt \
--GENOME_INDEX $REF \
--GENOME_LENGTH $gs \
--CHROM_SIZE_INFO $chrlength \
--CYTOBAND_DATA $cytoband \
--SPECIES 3 \
--output $OUT \
--prefix ${idx} \
--SELF_LIGATION_CUFOFF 8000 \
--EXTENSION_LENGTH 200 \
--MIN_COVERAGE_FOR_PEAK 5 \
--PEAK_MODE 1 \
--INPUT_ANCHOR_FILE ../01_seed_anchor/MH_H3K9ac_unify.IDR.narrowPeak \
--MIN_DISTANCE_BETWEEN_PEAK 500 \
--GENOME_COVERAGE_RATIO 0.8 \
--PVALUE_CUTOFF_PEAK 0.001 \
--PVALUE_CUTOFF_INTERACTION 0.05"
sleep 2
done