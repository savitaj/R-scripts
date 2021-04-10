# Immunarch is an R package for analysis of TCR/BCR repertoire data. 
# It allows you to look at shared clonotypes, clonal diversity, gene usage and more
# Supports all popular TCR and BCR analysis and post-analysis formats, including single-cell data: 
# ImmunoSEQ, IMGT, MiTCR, MiXCR, MiGEC, MigMap, VDJtools, tcR, AIRR, 10XGenomics, ArcherDX.

install.packages("devtools", dependencies = T)
devtools::install_github("immunomind/immunarch")

library("immunarch")
library("dtplyr")

data <-setwd("D:/Data files/Mixcr_Analysis_Partial_Assemble")
tcr<-repLoad(data)

#View(tcr$data)
names(tcr)
names(tcr$data)
top(tcr$data[1])

#Finding the number of shared clonotypes and visualize it as a correlation plot
ov = repOverlap(tcr$data, .method="public", .col = "aa", .a = 0.5, .b = 0.5, .verbose = T)
vis(ov)

# Using Shiny to make publication ready plots
p=vis(ov)
fixVis(p) 

# Visualize shared clonotypes as a heatmap with dendrogram
vis(ov, "heatmap2") 

# repOverlap function is designed to analyse the overlap between two or more repertoires
# .method = c("public", "overlap", "jaccard", "tversky", "cosine", "morisita", "inc+public", "inc+morisita")
# 'public' means shared clonotypes, "jaccard" index is percentage of objects two sets have in common out of total
imm_ov1 = repOverlap(tcr$data, .method = "public", .verbose = F)

#Morista is another statistical measure to check for the dispersion of individuals in a population. 
imm_ov2 = repOverlap(tcr$data, .method = "morisita", .verbose = F)
gridExtra::grid.arrange(vis(imm_ov1), vis(imm_ov2, .text.size=2), ncol = 2)

# Apply different analysis algorithms to the matrix of public clonotypes:
# "tsne" - t-Stochastic Neighbor Embedding
repOverlapAnalysis(imm_ov1, "tsne") %>% vis()

# Visualize CDR3 length distribution
repExplore(tcr$data, "lens") %>% vis()

# Repertoire Clonal Proportions
# Number of clonotypes occupying the 10% of repertoires
tcr_pr = repClonality(tcr$data, .method = "clonal.prop")
tcr_pr
vis(tcr_pr)

#Top clonal proportion
tcr_top = repClonality(tcr$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
tcr_top
vis(tcr_top)

#Rare clonal proportion
tcr_tail = repClonality(tcr$data, .method = "tail")
tcr_tail
vis(tcr_tail)
    
#Relative Abundance; summary of clonotypes with specific frequencies
tcr_hom = repClonality(tcr$data, .method = "homeo", .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)) 
tcr_hom
vis(tcr_hom) 

# Compare diversity of repertoires and visualise samples
# Chao1 estimator is a nonparameteric asymptotic estimator of species richness (number of species in a population)
div = repDiversity(tcr$data, .method = "chao1")
vis(div)

#visualize CDR3 length distribution
repExplore(tcr$data, "lens") %>% vis()

# The gene usage function allows to analyse the immune receptor gene usage 
# for (IGHD, IGHJ, IDHV, IGIJ, IGKJ, IGKV, IGLJ, IGLV, TRAJ, TRAV, TRBD, etc.) and gives the statistics
# geneUsage (
#  .data,
#  .gene = c("hs.trbv", "HomoSapiens.TRBJ", "macmul.IGHV"),
#  .quant = c(NA, "count"),
#  .ambig = c("inc", "exc", "maj"),
#  .type = c("segment", "allele", "family"),
#  .norm = F
#)


#Compute V gene usage and and highlight gene differences between different Age groups:
geneUsage(tcr$data[[1]]) %>% vis() 
#or
gu = geneUsage(tcr$data[[1]])
vis(gu, .by="Age", .meta=tcr$meta)
#TRBV gene usage
trbv_gu = geneUsage(tcr$data, "hs.trbv", .type = "family", .norm = T, .quant = "count", .ambig = "exc")
vis(trbv_gu, .plot = "hist", .grid = T)
#TRBD gene usage
trbd_gu = geneUsage(tcr$data, "hs.trbd", .type = "family", .norm = F, .quant = "count", .ambig = "exc")
vis(trbd_gu, .plot = "hist", .grid = T)
#TRAV gene usage
trav_gu = geneUsage(tcr$data, "hs.trav", .type = "family", .norm = T, .quant = "count", .ambig = "exc")
vis(trav_gu, .plot = "hist", .grid = T)
#TRAJ gene usage
traj_gu = geneUsage(tcr$data, "hs.traj", .type = "family", .norm = F, .quant = "count", .ambig = "exc")
vis(traj_gu, .plot = "hist", .grid = T)
gene_stats(traj_gu)

#cor - correlation coefficient for analysis. Can be "pearson", "kendall" or "spearman".
#Can be "pca", "mds", "js", "kmeans", "hclust", "dbscan" or "cor" if you want to calculate correlation coefficient.
#JS stands for Jensen-Shannon divergence, measuring the similarity between two probability distributions
tcr_gu_js = geneUsageAnalysis(trbv_gu, .method = "js", .verbose = F)
tcr_gu_cor = geneUsageAnalysis(trbv_gu, .method = "cor", .cor="pearson", .verbose = F)

gridExtra::grid.arrange(vis(tcr_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size=2.5, .signif.digits= 1), vis(tcr_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=2.5, .signif.digits= 2), ncol = 2)

# Track clonotypes
# Choose the first 10 amino acid clonotype sequences
# from the first repertoire to track
tc = trackClonotypes(tcr$data, list(1, 10), .col = "aa")
vis(tc, .plot = "smooth")

# Choose the first 10 clonotypes from the first repertoire
# with amino acid sequences and V segments
target = tcr$data[[1]] %>% select(CDR3.aa, V.name) %>% head(10)
tc = trackClonotypes(tcr$data, target)
vis(tc, .plot = "smooth")
vis(tc, .plot = "area")
vis(tc, .plot = "line")

pr.nt = pubRep(tcr$data, "nt", .verbose = F)
write.table(pr.nt, file="sharedclones.txt")

pr.aa = pubRep(tcr$data, "aa", .verbose = F)
write.table(pr.aa, file="Shared_clones.txt")

# Pass "aa+v" as the second parameter to build the public repertoire table using CDR3 aminoacid sequences and V alleles
# In order to use only CDR3 aminoacid sequences, just pass "aa"
pr.aav = pubRep(tcr$data, "aa+v", .verbose = F)
pr.aav

# You can also pass the ".coding" parameter to filter out all noncoding sequences first:
pr.aav.cod = pubRep(tcr$data, "aa+v", .coding=T)
pr.aav.cod

#explore differences in clonotype counts and volume

exp_cnt = repExplore(tcr$data, .method = "count")
vis(exp_cnt)

exp_vol = repExplore(tcr$data, .method = "volume")
vis(exp_vol)





