library(data.table)

#Expression quartiles

df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_deseq_sporo_v_other_output.txt")
df$quartile =   ntile(df$baseMean, 4)

#loop over expression quartiles and write site positions to file for other analyses (dfe alpha/regression etc)
for(i in c(1:4)) {
dfq = df[df$quartile == i,]


ids = fread("/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Genes_names_coordinates.txt", header = F)
ids[,5:7] = NULL


ids_quartile = ids[ids$V4 %in% dfq$row,]
ids_quartile$V4 = NULL

output_file = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/transcriptome/Marchantia_",i,"th_expression_quartile_sites.txt")
write.table(ids_quartile, , sep = "\t", col.names = F, row.names = F, quote = F)
 }



