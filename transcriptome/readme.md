# Files and scripts to perform gene expression analyses

## rna_pipeline.sh
Downloads raw transcriptomic data from JGI (requires login information to be added to jgi-query.config for JGI account!!!!)
Runs STAR to align data to genome, and HTSeq to get read counts for relevant genomic features (genes)

## R_run_Deseq.R
Takes input from HTSeq and runs DEseq 2. Produces normalized read counts/expression levels across all 4 tissues.

## expression_quartiles.R
Takes as input output from Deseq and calculates expression quartiles of genes based on normalized expression data, writes output for other scripts to use (regression analyses/DFE alpha etc)
