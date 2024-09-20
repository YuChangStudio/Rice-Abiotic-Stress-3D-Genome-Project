mkdir bam_QC
cd bam_QC
ln -s ../*rep*/*bam ./
ln -s ../*rep*/*bai ./

bsub -J bamsum -n 6 -R span[hosts=1]  -e e.e -q normal \
"multiBamSummary bins \
-b *.bam  \
-p 6 \
-bs 1000 \
--scalingFactors ChIApet.bam.scalingFactors \
--outRawCounts ChIApet.bam.count \
-o ./bam_summary.npz"