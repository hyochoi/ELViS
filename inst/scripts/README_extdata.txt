This file contains the explanation of the files in extdata

1. HPV16 Reference Files

The following reference files are included in extdata directory of ELViS.

inst/extdata/HPV16.fa
inst/extdata/HPV16.fa.fai
inst/extdata/HPV16REF_PaVE.gff

HPV16.fa is a reference sequence file with the associated fasta index file HPV16.fa.fai
HPV16REF_PaVE.gff are viral gene annotation file in GFF3 format.
They were downloaded from PaVE database (https://pave.niaid.nih.gov/locus_viewer?seq_id=HPV16REF)

2. BAM Files

The following BAM files are included in extdata directory of ELViS.

Control_100X_1102.bam
Control_100X_1102.bam.bai
Control_100X_1119.bam
Control_100X_1119.bam.bai

These are simulation data made using w-WESSIM-2. (https://github.com/GeorgetteTanner/w-Wessim2)
Briefly, sequencing reads were simulated using this tool and aligned to HPV16.fa to generate the BAM file.
