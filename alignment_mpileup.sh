bwamem2=/plas1/george.sandler/apps/bwaMEM2/bwa-mem2-2.0pre2_x64-linux/bwa-mem2

#prepare reference for bwa mem2
/plas1/george.sandler/apps/bwaMEM2/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index "/plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta"
java -jar /plas1/apps/picard/build/libs/picard.jar CreateSequenceDictionary R= /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta  O=/plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.dict


cd /plas1/george.sandler/marchantia_popgen_newgenome/

#Alignment
##############################


#ruderalif ref
parallel --jobs 8 '/plas1/george.sandler/apps/bwaMEM2/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta /plas2/george.sandler/marchantia_fastq/{}_R1.fastq.gz /plas2/george.sandler/marchantia_fastq/{}_R2.fastq.gz | samtools view -Sb -> /plas1/george.sandler/marchantia_popgen_newgenome/alignment/{}_bwa.bam' ::: Mar_D Mar_E Mar_G Mar_H Mar_J Mar_K Mar_L Mar_M Mar_N Mar_O Mar_S Mar_T Mar_U Mar_Y Mar_beta 




#Read processing
#############################
for ind in Mar_D Mar_E  Mar_G Mar_H Mar_J Mar_K Mar_L Mar_M Mar_N Mar_O Mar_S Mar_T Mar_U Mar_Y Mar_beta  
do

java -jar /plas1/apps/picard/build/libs/picard.jar AddOrReplaceReadGroups I= /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_bwa.bam O= /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_add.bam SO=coordinate RGID=${ind} RGLB=${ind} RGPL=Illumina RGPU=1 RGSM=${ind} TMP_DIR=/plas1/george.sandler/tmp/

java -jar /plas1/apps/picard/build/libs/picard.jar MarkDuplicates I=/plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_add.bam O= /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_ddpl.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M= /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_ddpl.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP= 3000 REMOVE_DUPLICATES=true

samtools index /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_ddpl.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M= /plas1/george.sandler/marchantia_popgen_newgenome/alignment/${ind}_ddpl.metrics

done


################################


#generate mpileup file for filtering
parallel --jobs 5 '/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools mpileup -d 1000 -q 50 -Q 20 -A -a AD,DP,ADF,ADR -Ov -f /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta /plas1/george.sandler/marchantia_popgen_newgenome/alignment/{}_ddpl.bam > /plas1/george.sandler/marchantia_popgen_newgenome/pileup/{}_bcftools.vcf' ::: Mar_D Mar_E Mar_G Mar_H Mar_J Mar_K Mar_L Mar_M Mar_N Mar_O Mar_S Mar_T Mar_U Mar_Y Mar_beta 
