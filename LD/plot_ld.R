

library(data.table)

#### LD script viz
library(data.table)


df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom_final.ld")
df = na.omit(df)
df = df[df$CHR_A == df$CHR_B,]

df$dist = df[,6]-df[,2]
df$dist2 = df$dist/1000
df$R2 = df$R*df$R

mod = smooth.spline(df$dist2,df$R2)


############## plot LD decay

plot(NULL, xlim=c(0,27000), ylim=c(0.12,0.38), ylab= expression(paste("r"^"2")), xlab="Distance between variants (kb)")
lines(mod, lwd=3)


