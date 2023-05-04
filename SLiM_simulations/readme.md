# files for running SLiM simulations
Basic workflow:

Edit SliM eidos file with desired parameters > run slim script in parallell > prepare output for hmmIBDD > run hmmIBD

## slim_admix_facsex.eidos slim_admix_facsex2k.eidos
Code for slim simulations with different rates of facultative sex (rate of sex editable within file).
First file for 10k pop sizes, second for 2k pop sizes

## slim_admix_lowrho.eidos
Code for slim simulations with different rates of recombination (rate recombination editable within file).

## run_slim.sh
bash file to run slim simulations in parallel

## prep_hmmibd.R
R script to take vcf output from SLiM and edit into input format for hmmIBD

## runhmmibd.sh
bash script to run hmmIBD on output of SLiM simulations in parallel
