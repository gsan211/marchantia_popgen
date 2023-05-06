#R inset to merge overlapping gene_windows
R
library(data.table)
library(GenomicRanges)

df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/genes_list_tabs_unique.list")

dt <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3")
r <- reduce(dt) #merge overlapping ranges
dc <- data.frame(r) #convert back to df
dc$width <- NULL
dc$strand <- NULL

write.table(dc, "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/genes_list_tabs_unique_collapsed.list", sep = "\t", quote = F, col.names =F, row.names =F)
q()
n
##################
#first tab seperated file of windows from R
cd  "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/"


bedtools intersect -a "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows" -b "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/genes_list_tabs_unique_collapsed.list" -wao > "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_coding_intersect.txt"

#calculate GC contect of windows
bedtools nuc -fi "/plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta" -bed "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows"   > "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_GC_content.txt"

############################################################################################################################################
#pull windows that overlap with sites of 4 expression quartiles and remove duplicate gene counts


bedtools intersect -a "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows" -b "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_1st_expression_quartile_sites.txt" -wo | sort -u > "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_quart1_intersect.txt"

############################################################################################################################################
#getting entropy/complexity statistics using Nessie


#first create multi fasta file with each window using bedtools
bedtools getfasta -fi /plas1/george.sandler/marchantia_popgen_newgenome/reference_new/Mpolylinkagemap_gapsestimated.fasta -bed /plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows > autosomes_10kb_windows.fa


#linguistic complexity
"/plas1/george.sandler/apps/nessie/nessie" -I "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/fastas/autosomes_10kb_windows.fa" -O /plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_lingcomplex.txt -L


#edite nessie output by moving chromosome id to same row as statistics and editing coordinates to be tab separated, also removes and records where N's in the reference cause 2 statistics to be reported
awk 'BEGIN{ ORS = "\t" }  /^>/ { print "\n", $0} NR>1{ print $2}' "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_lingcomplex.txt" \
| 
sed 's/>//g; s/-/\t/g; s/:/\t/g' | awk 'NF==4'> "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_lingcomplex_nomissing_edited.txt"
