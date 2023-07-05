# Data and scripts for running dfe alpha analysis

## prep_sfs_input.R
Takes as input data generated from SNPeff (see SNPeff scripts) and creates SFS for input to dfe_alpha
Have to add correct header to this file before inputting to dfe alpha, depends on SFS bins # etc. See dfe_alpha manual.


## run_dfe.sh
Runs dfe_alpha on SFS to generate NeS bins and estimates alpha/omega

## config files
Editable files which run_dfe.sh uses to feed commands into dfe_alpha. Need to be manually edited for each analysis. 

## SFS_unfolded_input / divergence_file
Input files used for dfe_alph analyses. Paths must be specified in config files. 

## bootstrap
Directory for bootstrapping SFS and running analysis to generate 95% confidence intervals for DFE_alpha estimates. 
Key file run_dfe_bootstrap.sh automatically bootstraps data and runs dfe_alpha on output
