
cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/newsubpop/"
parallel --jobs 20 /plas1/apps/SLiM.2.4.2/bin/slim -d o={} "/plas1/george.sandler/marchantia_popgen_newgenome/slim/slim_newsubpop.eidos"  ::: 1 2 3 4 5 



cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/"
parallel --jobs 20 /plas1/apps/SLiM.2.4.2/bin/slim -d o={} "/plas1/george.sandler/marchantia_popgen_newgenome/slim/slim_admix_subpop.eidos"  ::: 1 2 3 4 5 


cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/"
parallel --jobs 20 /plas1/apps/SLiM.2.4.2/bin/slim -d o={} "/plas1/george.sandler/marchantia_popgen_newgenome/slim/slim_admix_facsex.eidos" ::: 9 10 11 12 13 14 15


cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/"
parallel --jobs 20 /plas1/apps/SLiM.2.4.2/bin/slim -d o={} "/plas1/george.sandler/marchantia_popgen_newgenome/slim/slim_admix_facsex2k.eidos" ::: 2



cd "/plas1/george.sandler/marchantia_popgen_newgenome/slim/low_rho/"
parallel --jobs 20 /plas1/apps/SLiM.2.4.2/bin/slim -d o={} "/plas1/george.sandler/marchantia_popgen_newgenome/slim/slim_admix_lowrho.eidos" ::: 1 2 3 4 5 6 7 8

