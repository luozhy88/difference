
##############################lefser
library(lefser)
library(mia)
pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)


se <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)

se.rel <- lefser::relativeAb(se) 

res_group <- lefser(se.rel, groupCol = "group", kruskal.threshold = 0.05, lda.threshold = 0.0001, checkAbundances = FALSE)

lefserPlot(res_group, colors = c("red", "forestgreen"))
