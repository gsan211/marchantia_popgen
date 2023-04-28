
################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/"

parallel --jobs 7 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_facsex_90{9..15}

parallel --jobs 7 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o hmmibd/{}  ::: slim_facsex_90{9..15}



################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/newsubpop/"

parallel --jobs 5 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd_subsample.R" {} ::: slim_newsubpop_95_05_{1..5}

parallel --jobs 5 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o hmmibd/{}  ::: slim_newsubpop_95_05_{1..5}



################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/"

parallel --jobs 5 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_facsex2k_9999{1..5}

parallel --jobs 5 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o hmmibd/{}  ::: slim_facsex2k_00{1..5}


################################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/low_rho/"

parallel --jobs 8 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_rec_10_{1..8}

parallel --jobs 8 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o hmmibd/{}  ::: slim_rec_10_{1..8}





#################################
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/"

#parallel --jobs 7 Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" {} ::: slim_facsex_90{9..15}

parallel --jobs 7 "/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i  hmmibd/{}_input.txt -o test/{}  ::: slim_facsex2k_999{1..5}




Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd.R" slim_facsex2k_002

"/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i hmmibd/slim_facsex2k_902_input.txt  -o hmmibd/slim_facsex2k_902


cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/"
Rscript "/plas1/george.sandler/marchantia_popgen_newgenome/slim/prep_hmmibd_subsample.R" slim_facsex_906

"/plas1/george.sandler/apps/hmmibd/hmmIBD/hmmIBD_slim" -i hmmibd/slim_facsex_002_input.txt  -o hmmibd/slim_facsex_906
