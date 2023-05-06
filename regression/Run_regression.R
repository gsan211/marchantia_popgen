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

