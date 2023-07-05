library(data.table)
library(ggplot2)

df1 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_0fold.txt")
df1$cat = paste(df1$V2, df1$V4, sep = "_")
df1_sd = aggregate(df1, by=list(df1$cat), FUN= sd)
df1_sd$quart = "0 fold"

df2 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_3primeUTR.txt")
df2$cat = paste(df2$V2, df2$V4, sep = "_")
df2_sd = aggregate(df2, by=list(df2$cat), FUN= sd)
df2_sd$quart = "3 prime UTR"

df3 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_5primeUTR.txt")
df3$cat = paste(df3$V2, df3$V4, sep = "_")
df3_sd = aggregate(df3, by=list(df3$cat), FUN= sd)
df3_sd$quart = "5 prime UTR"

df4 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_intergenic.txt")
df4$cat = paste(df4$V2, df4$V4, sep = "_")
df4_sd = aggregate(df4, by=list(df4$cat), FUN= sd)
df4_sd$quart = "Intergenic"

df5 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_intron.txt")
df5$cat = paste(df5$V2, df5$V4, sep = "_")
df5_sd = aggregate(df5, by=list(df5$cat), FUN= sd)
df5_sd$quart = "Intronic"

df6 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_proportions_splice_site.txt")
df6$cat = paste(df6$V2, df6$V4, sep = "_")
df6_sd = aggregate(df6, by=list(df6$cat), FUN= sd)
df6_sd$quart = "Splice site"

boot = rbind(df1_sd, df2_sd, df3_sd, df4_sd, df5_sd, df6_sd)
boot$cat = boot$Group.1
boot$Group.1 = NULL
boot$V1 = NULL
boot$V2 = NULL
boot$V3 = NULL
boot$V4 = NULL
boot$V5 = NULL
boot$CI = 1.96*boot$V6

dat = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/True_DFE_props.txt", header =T)



dat2 <- merge(dat, boot, by = c("cat", "quart"))



ggplot(dat2, aes(x=cat, y=est,fill = as.factor(quart)))+ theme_bw() +theme(text = element_text(size = 18),  axis.text = element_text(colour = "black")) +
  geom_bar(position = 'dodge', stat='identity', width=.5,  colour="black")+  xlab("NeS category") + ylab("Proportion of sites")+
  geom_errorbar(aes(ymin = est+(CI), ymax = est-(CI)),
            width = 0.2,
            position = position_dodge(width = 0.5),
            color="black", size=0.5) +      
   scale_fill_manual(values=c("chocolate4", "chartreuse", "chartreuse4", "darkolivegreen", "darkgreen", "darkgoldenrod"), name=NULL) +   
   scale_x_discrete(labels=c("0_1" = "0-1", 
                             "1_10" = "1-10",
                             "10_100" = "10-100", 
                             "100_-99" = ">100")) +
    geom_hline(yintercept=0)                          

      
      
      
 ggsave(
  "DFE_plot_expanded_categories.pdf",
  plot = last_plot(),
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/",
  scale = 1,
  width = 20,
  height = 13,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
)           
            
            
 

