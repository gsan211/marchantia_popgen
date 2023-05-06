library(modelsummary)
library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(gridExtra)

win6 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/regression_input.txt", header = T)


#model for pi
modpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) +  scale(ling_complex) + scale(exp) + scale(pi_n) , data = win6)
summary(modpi)
modelplot(modpi)

#model for d_xy
moddxy <- lm(d_xy_meanm ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) + scale(ling_complex) + scale(exp) + scale(d_xy_n), data = win6)
summary(moddxy)
modelplot(moddxy)


models = list(modpi, moddxy)

modelplot(models)



##########################
#models to explore how callable site# affects regression results
#fit separate models for low vs high callabel sites windows
lowcal = win6[win6$pi_n < 1200,]
modlowpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds)  + scale(ling_complex) + scale(exp) + scale(pi_n) + scale(frac_GC) , data = lowcal)

highcal = win6[win6$pi_n > 2100,]
modhighpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(ling_complex) + scale(exp) + scale(pi_n) + scale(frac_GC) , data = highcal)


models = list(modlowpi, modhighpi)

modelplot(models)

#########################################################################
#bootstrapping p values
library(boot.pval)


win7 = win6
win7$d_xy_meanm = NULL
win7$d_xy_n = NULL
win7$win_mid = NULL
win7$entropy = NULL
win7$id = NULL


win7 = na.omit(win7)
modpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) +  scale(ling_complex) + scale(exp) + scale(pi_n) , data = win7)
#modpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) +  scale(ling_complex) + scale(exp) + scale(pi_n) + scale(entropy) , data = win6)
summary(modpi)


boot_summary(modpi, R = 999)

########

moddxy <- lm(d_xy_meanm ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) + scale(ling_complex) + scale(exp) + scale(d_xy_n), data = win6)
summary(moddxy, R = 999)

win9 = win6
win9$winpi = NULL
win9$pi_n = NULL
win9$win_mid = NULL
win9$entropy = NULL
win9$id = NULL
win9 = na.omit(win9)

moddxy <- lm(d_xy_meanm ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(frac_GC) + scale(ling_complex) + scale(exp) + scale(d_xy_n), data = win9)
summary(moddxy)

boot_summary(moddxy)

############################################


lowcal = win7[win7$pi_n < 1200,]
modlowpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds)  + scale(ling_complex) + scale(exp) + scale(pi_n) + scale(frac_GC) , data = lowcal)

highcal = win7[win7$pi_n > 2100,]
modhighpi <- lm(winpi ~ V1 + scale(frac_utr) + scale(frac_cds) + scale(ling_complex) + scale(exp) + scale(pi_n) + scale(frac_GC) , data = highcal)

summary(modlowpi)

summary(modhighpi)


boot_summary(modlowpi, R = 999)


boot_summary(modhighpi, R = 999)


