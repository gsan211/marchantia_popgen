
################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/"

parallel --jobs 7 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_facsex_90{1..5}

parallel --jobs 7 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o hmmibd/{}  ::: slim_facsex_90{1..5}




#################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/"

parallel --jobs 7 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_facsex_90{1..5}

parallel --jobs 7 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o test/{}  ::: slim_facsex2k_999{1..5}


