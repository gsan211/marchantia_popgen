cd "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/"

#download files from JGI, need to enter own credential in jgi-query.config!!!!!!!!!!!!
python "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/jgi-query.py " Marpolnscriptome

1:1-8



#split files

cd /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/fastq/
for ind in sporophyte thallus archegonia antheridia
do

./deinterleave_fastq.sh < Marchantia_${ind}_RNA.fastq  Marchantia_${ind}_RNA_R1.fastq  Marchantia_${ind}_RNA_R2.fastq
done



#generate index for STAR
#/plas1/george.sandler/apps/STAR-2.6.0a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/ --genomeFastaFiles /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta  --runThreadN 15

#run star
for ind in Marchantia_antheridia Marchantia_archegonia Marchantia_sporophyte Marchantia_thallus
do
/plas1/george.sandler/apps/STAR-2.6.0a/bin/Linux_x86_64/STAR --runThreadN 15 --genomeDir /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/ --readFilesIn /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/fastq/${ind}_RNA_R1.fastq /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/fastq/${ind}_RNA_R2.fastq --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/star_aligned/${ind}

done


#HTSeq

htseq-count --nonunique random --type=CDS --idattr=Parent --additional-attr=Name /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/star_aligned/Marchantia_archegoniaAligned.sortedByCoord.out.bam /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.gff3 > /plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_archegonia_htseq_counts.txt


#scrape gene names from gff file, could modify to retain coordinates too

#awk {$1 =$1}1 means for non-zero records set them to same value, necessary otherwise awk won't edit output field seperator as requested. The extra "1" specifies non-zero values otherwise awk will chuck them (?)
awk '$3 == "mRNA"' /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated_clean.gff3 | awk '{print $1,$4,$5,$9}'| awk 'BEGIN {OFS=";"} { $1=$1 }1'|   awk -F';' 'BEGIN {OFS="\t"} { $1=$1 } 1'  | sed 's/ID=//g' |sed 's/Name=//g' > "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Genes_names_coordinates.txt"


#paste together the HTSeq output for each tissue
cd "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/"
paste Marchantia_antheridia_htseq_counts.txt Marchantia_archegonia_htseq_counts.txt Marchantia_sporophyte_htseq_counts.txt Marchantia_thallus_htseq_counts.txt >  Marchantia_htseq_counts.txt 















