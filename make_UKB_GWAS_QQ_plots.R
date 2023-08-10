# set working directory for dev and debug
bdir = "/Users/williamlee/SharpLab/"
setwd(bdir)

library(ggplot2)

###############################################################################

# post-meta-analysis GWAS

gwas_df = read.table(sep="\t",header=T,file='out/UKB_GWAS_meta1.txt', stringsAsFactors=FALSE)

pvector = gwas_df$P.value
o = -log10(sort(pvector, decreasing = FALSE))
e = -log10(ppoints(length(pvector)))
p = data.frame(observed = o, expected = e)

g = ggplot(p, aes(expected, observed)) +
  geom_abline(intercept = 0, slope = 1, size = 0.4) +
  geom_point() +
  xlab("Expected P value (-log10)") +
  ylab("Observed P value (-log10)") +
  theme_bw(base_size = 20)

jpeg("out/UKB_GWAS_meta_pval_QQ_plot.jpeg", width = 2400, height = 2400, res = 300, pointsize = 12)
print(g)
dev.off()

###############################################################################

# SC GWAS

gwas_df = read.table(sep=" ",header=T,file='out/UKB_GWAS_SC_pval.txt', stringsAsFactors=FALSE)

gwas_df$P.value = 10**-(gwas_df$LOG10P)

pvector = gwas_df$P.value
o = -log10(sort(pvector, decreasing = FALSE))
e = -log10(ppoints(length(pvector)))
p = data.frame(observed = o, expected = e)

g = ggplot(p, aes(expected, observed)) +
  geom_abline(intercept = 0, slope = 1, size = 0.4) +
  geom_point() +
  xlab("Expected P value (-log10)") +
  ylab("Observed P value (-log10)") +
  theme_bw(base_size = 20)

jpeg("out/UKB_GWAS_SC_pval_QQ_plot.jpeg", width = 2400, height = 2400, res = 300, pointsize = 12)
print(g)
dev.off()

###############################################################################

# deCODE GWAS

gwas_df = read.table(sep=" ",header=T,file='out/UKB_GWAS_deCODE_pval.txt', stringsAsFactors=FALSE)

gwas_df$P.value = 10**-(gwas_df$LOG10P)

pvector = gwas_df$P.value
o = -log10(sort(pvector, decreasing = FALSE))
e = -log10(ppoints(length(pvector)))
p = data.frame(observed = o, expected = e)

g = ggplot(p, aes(expected, observed)) +
  geom_abline(intercept = 0, slope = 1, size = 0.4) +
  geom_point() +
  xlab("Expected P value (-log10)") +
  ylab("Observed P value (-log10)") +
  theme_bw(base_size = 20)

jpeg("out/UKB_GWAS_deCODE_pval_QQ_plot.jpeg", width = 2400, height = 2400, res = 300, pointsize = 12)
print(g)
dev.off()

###############################################################################