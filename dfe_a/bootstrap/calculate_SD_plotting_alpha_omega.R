library(data.table)
library(ggplot2)

df1 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_0fold.txt")
df1_sd = aggregate(df1, by=list(df1$V1), FUN= sd)
df1_sd$quart = "0 fold"

df2 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_3primeUTR.txt")
df2_sd = aggregate(df2, by=list(df2$V1), FUN= sd)
df2_sd$quart = "3 prime UTR"

df3 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_5primeUTR.txt")
df3_sd = aggregate(df3, by=list(df3$V1), FUN= sd)
df3_sd$quart = "5 prime UTR"

df4 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_intergenic.txt")
df4_sd = aggregate(df4, by=list(df4$V1), FUN= sd)
df4_sd$quart = "Intergenic"

df5 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_intron.txt")
df5_sd = aggregate(df5, by=list(df5$V1), FUN= sd)
df5_sd$quart = "Intronic"

df6 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/DFE_bootstrapped_alpha_omega_splice_site.txt")
df6_sd = aggregate(df6, by=list(df6$V1), FUN= sd)
df6_sd$quart = "Splice site"

boot = rbind(df1_sd, df2_sd, df3_sd, df4_sd, df5_sd, df6_sd)
boot$Group.1 = NULL
boot$V1 = NULL
boot$V3 = NULL
boot$V5 = NULL
boot$CI_adapt_diverence = 1.96*boot$V2
boot$CI_alpha = 1.96*boot$V4
boot$CI_omega = 1.96*boot$V6

dat = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/True_omega_alpha.txt", header =T)



dat2 <- merge(dat, boot, by = c("quart"))



ggplot(dat2, aes(x=quart, y=alpha))+ theme_bw() +geom_point(position = 'dodge', stat='identity', width=.5)+ 
  geom_errorbar(aes(ymin = alpha+(CI_alpha), ymax = alpha-(CI_alpha)),
            width = 0.2,
            position = position_dodge(width = 0.5),
            color="red", size=0.5)
      

ggplot(dat2, aes(x=quart, y=omega))+ theme_bw() +geom_point(position = 'dodge', stat='identity', width=.5)+ 
  geom_errorbar(aes(ymin = omega+(CI_omega), ymax = omega-(CI_omega)),
            width = 0.2,
            position = position_dodge(width = 0.5),
            color="red", size=0.5)
      

#######################
dat3 = dat2[dat2$quart == "0 fold" | dat2$quart == "3 prime UTR"| dat2$quart == "5 prime UTR",]   


fin = ggplot(dat3, aes(x=quart, y=alpha, colour = quart))+ theme_bw() + geom_point(colour = "black", size = 3.9) + geom_point(position = 'dodge', stat='identity', width=.5, size=3)+ theme(text = element_text(size = 18),axis.text=element_text(colour = "black"), panel.grid.minor = element_blank() , panel.grid.major = element_blank())+ xlab("Site Type") + ylab("Alpha") + ylim(-0.1, 0.6)+ scale_fill_manual(values=c("darkgreen","darkred","red")) + 
  geom_errorbar(aes(ymin = alpha+(CI_alpha), ymax = alpha-(CI_alpha)),
            width = 0.2,
            position = position_dodge(width = 0.5),
            color="black", size=0.5) + 
        scale_color_manual(values = c("0 fold" = "coral4",
                                "3 prime UTR"="chartreuse",
                                "5 prime UTR"="darkolivegreen3")) + ylim(-0.1, 0.6)


pl1 = ggplot(df_all, aes(x=fract_sites_IBD, y = ..scaled.., fill=sbsp)) + theme_bw() + geom_density(alpha=0.3) +
    theme(text = element_text(size = 18),axis.text=element_text(colour = "black"),legend.position = c(0.7, 0.7), panel.grid.minor = element_blank() , panel.grid.major = element_blank()) +
   xlab("Fraction of genome IBD per sample pair") + ylab("Scaled Density") +
   scale_x_continuous(limits = c(0, 0.6))+
   scale_fill_manual(values=c("darkgreen","darkred"), name="Subspecies",  labels = c(expression(italic("ruderalis")), expression(italic("polymorpha"))  )) 
 

ggsave(
  "alpha_SD.pdf",
  plot = fin,
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/",
  scale = 1,
  width = 13,
  height = 13,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
)           
            
           