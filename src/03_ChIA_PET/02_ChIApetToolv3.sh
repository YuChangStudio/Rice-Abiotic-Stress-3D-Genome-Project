for i in *.gz
do
bsub -n 1 -q normal \
"
gzip -d $i
"
sleep 2
done

REF=~/genome_anno/mh63/bwaid/MH63RS3bwaid
chrlength=~/genome_anno/mh63/MH63RS3.chrlength.txt
cytoband=~/genome_anno/mh63/MH63RS3.cytoBandIdeo.txt
gs=3.9E8
OUT=../03_out_spcific_region/

for i in H3K9ac_D H3K9ac_N H3K9ac_RE
do
for j in rep1 rep2 pooled
do
bsub -n 8 -q q2680v2 -o cptv3.o -e cptv3.e -R "rusage[mem=30GB]" -R span[hosts=1] \
"java -jar ~/software/ChIA-PET_Tool_V3-master/ChIA-PET.jar \
--mode 1 \
--thread 8 \
--start_step 1 \
--fastq1 MH_${i}_ChIApet_${j}_1.clean.fq \
--fastq2 MH_${i}_ChIApet_${j}_2.clean.fq \
--minimum_linker_alignment_score 19 \
--linker ~/software/ChIA-PET_Tool_V3-master/linker/linker_long.txt \
--GENOME_INDEX $REF \
--GENOME_LENGTH $gs \
--CHROM_SIZE_INFO $chrlength \
--CYTOBAND_DATA $cytoband \
--SPECIES 3 \
--output $OUT \
--prefix MH_${i}_${j} \
--SELF_LIGATION_CUFOFF 8000 \
--EXTENSION_LENGTH 500 \
--MIN_COVERAGE_FOR_PEAK 5 \
--PEAK_MODE 1 \
--INPUT_ANCHOR_FILE ../02_spanPeak/MH_${i}.span.narrowPeak \
--MIN_DISTANCE_BETWEEN_PEAK 500 \
--GENOME_COVERAGE_RATIO 0.8 \
--PVALUE_CUTOFF_PEAK 0.001 \
--PVALUE_CUTOFF_INTERACTION 0.05"
sleep 2
done
done

for i in MH*
do
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/${i}/${i}.ipet ./${i}/${i}.ipet
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/${i}/${i}.spet ./${i}/${i}.spet
ln -s ~/3D_project/02_MH_ChIApet/03_out_spcific_region/${i}/${i}.opet ./${i}/${i}.opet
done

REF=~/genome_anno/mh63/bwaid/MH63RS3bwaid
chrlength=~/genome_anno/mh63/MH63RS3.chrlength.txt
cytoband=~/genome_anno/mh63/MH63RS3.cytoBandIdeo.txt
gs=3.9E8
OUT2=../04_out_unify_region/01_ChIAPETtoolv3_out

for i in H3K9ac_D H3K9ac_N H3K9ac_RE
do
for j in rep1 rep2 pooled
do
bsub -n 1 -q q2680v2 -o cptv3.o -e cptv3.e -R "rusage[mem=30GB]" -R span[hosts=1] \
"java -jar ~/software/ChIA-PET_Tool_V3-master/ChIA-PET.jar \
--mode 1 \
--thread 8 \
--start_step 5 \
--fastq1 MH_${i}_ChIApet_${j}_1.clean.fq \
--fastq2 MH_${i}_ChIApet_${j}_2.clean.fq \
--minimum_linker_alignment_score 19 \
--linker ~/software/ChIA-PET_Tool_V3-master/linker/linker_long.txt \
--GENOME_INDEX $REF \
--GENOME_LENGTH $gs \
--CHROM_SIZE_INFO $chrlength \
--CYTOBAND_DATA $cytoband \
--SPECIES 3 \
--output $OUT2 \
--prefix MH_${i}_${j} \
--SELF_LIGATION_CUFOFF 8000 \
--EXTENSION_LENGTH 500 \
--MIN_COVERAGE_FOR_PEAK 5 \
--PEAK_MODE 1 \
--INPUT_ANCHOR_FILE ../02_unified_peaks/MH_H3K9ac_unify.IDR.narrowPeak \
--MIN_DISTANCE_BETWEEN_PEAK 500 \
--GENOME_COVERAGE_RATIO 0.8 \
--PVALUE_CUTOFF_PEAK 0.001 \
--PVALUE_CUTOFF_INTERACTION 0.05"
sleep 2
done
done
