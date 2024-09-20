bsub -n 1 -e e.e -q q2680v2 \
"pyGenomeTracks \
--tracks pyGenometrack.Fig5n.ini \
--region Chr11:17,980,000-18,020,000 \
--dpi 600 \
-o Fig5n.track.pdf"

bsub -n 1 -e e.e -q q2680v2 \
"pyGenomeTracks \
--tracks pyGenometrack.Fig2c.ini \
--region Chr03:27300000-27685000 \
--dpi 600 \
-o Fig2c.track.pdf"