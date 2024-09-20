#detect motifs in DS specific SPRs
ame --verbose 1 --oc . --scoring totalhits --method fisher --hit-lo-fraction 0.1 --evalue-report-threshold 100.0 --control N_dominate_superAnchor.span500.fa D_dominate_superAnchor.span500.fa motif_db/JASPAR/JASPAR2022_CORE_plants_non-redundant_v2.meme

#the bzip motif enrichment in the bzip23 binding sites.
centrimo --oc . --verbosity 1 --score 5.0 --ethresh 10.0 --bfile bzip23.DS.fa.bg bzip23.DS.fa versions.meme