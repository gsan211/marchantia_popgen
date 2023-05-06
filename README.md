# marchantia_popgen
# Scripts associated with analysis of M. polymorpha population genomic data

# Core pipeline scripts

## alignment_mpileup.sh
Script for aligning data, read processing and generating bcftools mpileup files for each sample

## individual_filtering
Scripts for filtering of mpileup files and generating combined VCF of genotype calls

## R_calculate_pi_dxy.R
Script for calculating pi based on genotype table generated after all filtering done in previous steo

## R_calculate_pi_within _sbsp.R
Scripts for calculating d_xy between subspecies



#####################


## Downstream analyses
## SLiM_simualtions
Scripts for running simulatotions and calculating IBD statistics

## regression
Scripts and file inputs for running regression of various genomic featurs and diversity/diverence
