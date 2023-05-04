# Master files for filetring and calling SNPs used in dataset

Files setup for filtering ssp ruderalis, but same workflow used for sbsp polymorpha

## Main workflow
Main workflow is outline in file:
"/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/R_processing_prepy.R"

Briefly
1) BCFtools pileup file is used to extract allelic read counts at each site for each individual, R script used to split those counts up into separate columns
2) Python script run to determine # problematic sites in each genomic window for each individual
3) R script run to pinpoint especially bad windows for each individual, and concatenate them together, output written
4) BCFtools used to generate genotype calls for each individual
5) Filtering then applied to those calls for each individual (including removing bad windows for that individual)
6) All individuals' genotype calls merged, and group level filtering applied

