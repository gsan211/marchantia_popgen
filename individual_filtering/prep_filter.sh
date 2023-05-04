#      "/plas1/george.sandler/apps/htslib/bgzip" 
#"/plas1/george.sandler/apps/htslib/tabix"
bcftoolsg=/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools
samtoolsg=/plas1/george.sandler/apps/samtools-1.8/samtools



#parse allelic depth information from pileup file for window analysis
parallel --jobs 15 '/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools query /plas1/george.sandler/marchantia_popgen_newgenome/pileup/{}_bcftools.vcf -f "%CHROM\t%POS\t%ALT\t[\t%AD]\n" -o /plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/Mar_bcftools_AD_{}.txt'  ::: Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_E Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O 

#run R script to separate out allelic depths into different columns from BCFtools output

parallel --jobs 15 'Rscript /plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/R_processing_prepy.R {}' ::: Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_E Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O 

#####################################################################################
#individual window based filtering


#running the python script to first calculate # bad sites in windows, and then filter those windows using R script 2
#python script run on each chromosome separately, chromosomes name taken as arg2 

parallel --jobs 15 'python /plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/python/Py_genomic_window_filtering.py {} chr1'  :::  Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O Mar_E Mar_U



#merge separate windows output together (merges chromosomes)
#cd "/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/"
for ind in Mar_E Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O
do
cat /plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/${ind}_*_windows.txt >  /plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/${ind}_allwindows.txt

done 

#takes in window input, chooses bad windows (2 or more bad sites in individual), merges all bad windows together and writes output file for filtering in BCFtools
parallel --jobs 7 '/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/R_processing2.R {}'  :::  Mar_E Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O 




####################################################################################


#call SNPs 

parallel --jobs 15 '/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools call -Ov -m /plas1/george.sandler/marchantia_popgen_newgenome/pileup/{}_bcftools.vcf --ploidy 1  > /plas1/george.sandler/marchantia_popgen_newgenome/calls/{}_calls.vcf'  ::: Mar_E Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O 


#do individual level filtering on SNPs including coverage, and removing bad regions 
#remove indels, low coverage sites within samples, bad regions
for ind in  Mar_E Mar_U Mar_S Mar_J Mar_N Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_G Mar_H Mar_K Mar_M Mar_O
do

/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools filter -Ov /plas1/george.sandler/marchantia_popgen_newgenome/calls/${ind}_calls.vcf -e 'TYPE = "indel"| INFO/DP<3' -T ^/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/${ind}_bad_windows.txt> /plas1/george.sandler/marchantia_popgen_newgenome/calls/${ind}_calls_filtered.vcf

/plas1/george.sandler/apps/htslib/bgzip  /plas1/george.sandler/marchantia_popgen_newgenome/calls/${ind}_calls_filtered.vcf
/plas1/george.sandler/apps/htslib/tabix /plas1/george.sandler/marchantia_popgen_newgenome/calls/${ind}_calls_filtered.vcf.gz


done 





#merge individuals
cd "/plas1/george.sandler/marchantia_popgen_newgenome/calls/"
/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools merge -Ov /plas1/george.sandler/marchantia_popgen_newgenome/calls/*_calls_filtered.vcf.gz > /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered.vcf


#subspecies level filtering
#filter sites where more than 2/3 samples is missing, and sites where total depth over 80 polymorpha

$bcftoolsg filter -Ov -e 'TYPE = "indel"' /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered.vcf|\
$bcftoolsg filter -Ov -e 'N_ALT >1' | \
$bcftoolsg filter -Ov -e 'F_MISSING > 0.32' | \
$bcftoolsg filter -Ov -e 'F_PASS(FORMAT/DP< 1) > 0.32' | \
$bcftoolsg filter -Ov -e 'INFO/DP < 5' | \
$bcftoolsg filter -Ov -e 'INFO/DP > 200' >  /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.vcf

#query output genotypes to usable table for downstream diversity estimates etc.
$bcftoolsg query /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.vcf -f "%CHROM\t%POS\t%ALT\t[\t%GT]\n" -o /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.GTS.txt












######################################################
######################################################
#sex chromosome specific filtering

#males
#parallel --jobs 8 '/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools filter -Ov -t Chr_Y_A,Chr_Y_B /plas1/george.sandler/marchantia_popgen_newgenome/calls/{}_calls_filtered.vcf.gz > /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/{}_V_filtered.vcf' :::  Mar_E Mar_U Mar_S Mar_J Mar_N Mar_E Mar_G Mar_M 




#for ind in Mar_E Mar_U Mar_S Mar_J Mar_N Mar_E Mar_G Mar_M 
#do

#/plas1/george.sandler/apps/htslib/bgzip  /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/${ind}_V_filtered.vcf
#/plas1/george.sandler/apps/htslib/tabix  /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/${ind}_V_filtered.vcf.gz

#done

#/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools merge -Ov /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/*_V_filtered.vcf.gz > /plas1/george.sandler#/marchantia_popgen_newgenome/calls/UV/Mar_all_V_filtered.vcf &


###############################################

#females
parallel --jobs 8 '/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools filter -Ov -t scaffold_17,scaffold_18,scaffold_210,scaffold_227,scaffold_230,scaffold_240,scaffold_250,scaffold_277,scaffold_497 /plas1/george.sandler/marchantia_popgen_newgenome/calls/{}_calls_filtered.vcf.gz > /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/{}_U_filtered.vcf' :::   Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_H Mar_K Mar_O  


#for ind in Mar_D Mar_beta Mar_T Mar_L Mar_Y Mar_H Mar_K Mar_O 
#do

#/plas1/george.sandler/apps/htslib/bgzip  /plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/filtered_bcftools_UV/${ind}_U_filtered.vcf
#/plas1/george.sandler/apps/htslib/tabix  /plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/filtered_bcftools_UV/${ind}_U_filtered.vcf.gz

#done

#/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools merge -Ov /plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/filtered_bcftools_UV/*_U_filtered.vcf.gz > /plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/filtered_bcftools_UV/Mar_polymorpha_U_filtered.vcf.gz


###############################
###############################
$bcftoolsg filter -Ov -e 'TYPE = "indel"' /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_V_filtered.vcf |\
$bcftoolsg filter -Ov -e 'N_ALT >1' | \
$bcftoolsg filter -Ov -e 'F_MISSING > 0.32' | \
$bcftoolsg filter -Ov -e 'F_PASS(FORMAT/DP< 1) > 0.32' | \
$bcftoolsg filter -Ov -e 'INFO/DP < 4' | \
$bcftoolsg filter -Ov -e 'INFO/DP > 120' > /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_V_secondfiltered.vcf

$bcftoolsg query /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_V_secondfiltered.vcf -f "%CHROM\t%POS\t%ALT\t[\t%GT]\n" -o /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_V_secondfiltered.GTS.txt






