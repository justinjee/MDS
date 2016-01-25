#
MDS

Code is organized as follows:
1) parsefastq.py takes a fastq file (from an MDS experiment) and breaks it up into several smaller files by ROI and phase
2) fastq2snps.py and fastq2frameshifts.py (somewhat improperly named) takes each of these smaller files and produces a list of variants

Data is organized as follows:
Most recent run, 20151215, contains data on rpoB reverse complement, mrcA (3'), dual-barcoded exogenous sequence, and mrcA in high translation strain and high transcription strain
All other files are annotated as in the SRA repo (SRP64302) 

NOTE: folders are organized according to number of reads, where paired-end reads count as 2 reads. So strict-3, for example, corresponds to data in which R (the number of separate DNA molecules with the same barcode) = 2.
