
#####################################################
#####################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/"



for site in 3primeUTR 5primeUTR 0fold splice_site intron intergenic
do

for i in {1..1000}
do

#r script takes original sites and bootstraps them into a new scrambled SFS
Rscript /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/bootstrap_data.R $site

#add correct header to SFS for DFE_a input
cat "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/header_sfs" "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/bootstrapped_sfs" > "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/bootstrapped_sfs_input"

#first run on neutral SFS to estimate demography
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/est_dfe -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/config_file_epoch2_neut.txt

#then on selected
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/est_dfe -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/config_file_epoch2_sel.txt

#calculate Nes proportions and save to file for each site type
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/prop_muts_in_s_ranges -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/sel/est_dfe.out -o dfe_output_bootstrap | grep -e "lower" >> DFE_bootstrapped_proportions_${site}.txt



#estimate alpha/omega

/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/est_alpha_omega -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/config_est_alpha_omega.txt| grep -e 'adaptive_divergence' >>  DFE_bootstrapped_alpha_omega_${site}.txt



done
done

