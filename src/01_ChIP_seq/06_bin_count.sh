bsub -J bamsum -n 6 -R span[hosts=1] -e e.e -q q2680v2 \
"multiBamSummary bins \
-b *rmdup.bam \
-p 6 \
-bs 10000 \
--scalingFactors ./QC/bam.10k.scalingFactors \
--outRawCounts ./QC/bam.10k.count \
-o ./QC/bam_summary.10k.npz"