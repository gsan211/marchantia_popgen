
#####################################################
#####################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/"

#first run on neutral SFS to estimate demography
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/est_dfe -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/config_file_epoch2_neut.txt

#then on selected
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/est_dfe -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/config_file_epoch2_sel.txt

#
/plas1/george.sandler/apps/dfe_alpha/dfe-alpha-release-2.16/prop_muts_in_s_ranges -c /plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/sel/est_dfe.out -o CHECK
