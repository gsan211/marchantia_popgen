library(data.table)
library(DESeq2)

df = read.table("/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_htseq_counts.txt", row.names=1, sep = "\t", header = T) #ignore warning
df = df[1:(nrow(df)-5),] #remove last 5 rows which have non-genic counts (non genic reads, poor quality reads etc)


df$X = NULL
df$X.1 = NULL
df$X.2 = NULL
df$X.3 = NULL

df$archegonia = NULL
df$sporophyte = NULL
df$thallus = NULL

#df$archegonia_count = NULL
#df$antheridia_count = NULL

#set treatments here, run once for every tissue as the focal group in analysis, sporophte shown as exmaple
tissue <- factor(c(rep("other",2),rep("sporo", 1), rep("other",1)))    #Set treatments

coldata <- data.frame(names(df)[-5], tissue)

dds <- DESeqDataSetFromMatrix(countData=df, colData=coldata, design=~tissue, tidy=FALSE)
#dds

dds <- DESeq(dds)             #RUN DESEQ2

#negative log fold = higher sporo exp


res <- results(dds, tidy=TRUE)
#head(results(dds, tidy=TRUE))

#output results from all sites for expression quartiles analysis
write.table(ids_sporo, "/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_deseq_sporo_v_other_output.txt", sep = "\t", col.names = F, row.names = F, quote = F)

res <- res[order(res$padj),]

res2 = res[res$padj <0.05,]
res2 = na.omit(res2)


sporo_over = res2[res2$log2FoldChange > 0 ,]

###########################
ids = fread("/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Genes_names_coordinates.txt", header = F)
ids[,5:7] = NULL



ids_sporo = ids[ids$V4 %in% sporo_over$row,]
ids_sporo$V4 = NULL

#significantly differntialy expressed genes
write.table(ids_sporo, "/plas1/george.sandler/marchantia_popgen_newgenome/revisions/dfe_sporo/Marchantia_sporo_overexpressed_FDR005_sites.txt", sep = "\t", col.names = F, row.names = F, quote = F)

