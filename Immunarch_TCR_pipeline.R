install.packages("devtools", dependencies = T)
install.packages("immunarch")
install.packages("tcR")


library("immunarch")
library("tcR")

setwd("~/MIXCR_pipeline1/")
tcr <-repLoad("~/MIXCR_Pipeline1/")

#View(tcr$data)
names(tcr)
names(tcr$data)
top(tcr$data[1])

ov = repOverlap(tcr$data, .method="public", .col = "aa", .quant = "count", .a = 0.5, .b = 0.5, .verbose = T, .dup = "merge")
vis(ov, "heatmap2")


imm_ov1 = repOverlap(tcr$data, .method = "public", .verbose = F)
imm_ov2 = repOverlap(tcr$data, .method = "morisita", .verbose = F)
grid.arrange(vis(imm_ov1), vis(imm_ov2, .text.size=2), ncol = 2)

vis(imm_ov1, "heatmap2")

# Multi-dimensional Scaling
vis(repOverlapAnalysis(imm_ov1, "mds"))

#t-stocastic neighbour embedding
vis(repOverlapAnalysis(imm_ov1, "tsne"))

# Clusterise the MDS resulting components using K-means
vis(repOverlapAnalysis(imm_ov1, "mds+kmeans"))

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

tcr_pr = repClonality(tcr$data, .method = "clonal.prop")
vis(tcr_pr)
p=vis(tcr_pr)
fixVis(p)
vis(tcr_pr, .by=c("parameter1", "parameter2"), .meta=tcr$meta)
p=vis(tcr_pr, .by=c("RNA","pipeline"), .meta=tcr$meta)
fixVis(p)

tcr_top = repClonality(tcr$data, .method = "top", .head = c(10, 100, 1000, 3000))
vis(tcr_top)
#vis(tcr_top, .by=c("pipeline", "RNA"), .meta=tcr$meta)
p=vis(tcr_top)
fixVis(p)

tcr_tail = repClonality(tcr$data, .method = "tail")
tcr_tail

#Clonotype abundance plot using Shiny
tcr_hom = repClonality(tcr$data, .method = "homeo", .clone.types = c(Small = .0001, Medium = .001, Large = .01, Expanded = 0.1, Hyperexpanded = 1)) 
vis(tcr_hom)
vis(tcr_hom, .by=c("pipeline", "RNA"), .meta=tcr$meta)
hom = vis(tcr_hom)
fixVis(hom)

#vis(tcr_top) #plots the summary proportion of clonotypes with specific indices
#vis(tcr_tail) #tail clonal proportion
#vis(tcr_hom)

#explore differences in clonotype length
exp_len = repExplore(tcr$data, .method = "len", .col = "aa")
vis(exp_len)

exp_cnt = repExplore(tcr$data, .method = "count")
vis(exp_cnt)

exp_vol = repExplore(tcr$data, .method = "volume")
vis(exp_vol)


tcr_gu = geneUsage(tcr$data)


gu.clust = geneUsageAnalysis(tcr_gu, .method = "js+hclust")
vis(gu.clust)
gu.kmeans = geneUsageAnalysis(tcr_gu, .method = "pca+kmeans")
vis(gu.kmeans)

gu.kmeans = geneUsageAnalysis(tcr_gu, .method = "pca+kmeans", .base = 2, .norm.entropy = F, .cor = "pearson", .do.norm = T, .laplace = 1e-12, .verbose = T, .k = 2, .eps = 0.01, .perp = 1, .theta = 0.1)
vis(gu.kmeans)
gu.kmeans = geneUsageAnalysis(tcr_gu, .method = "js+pca+kmeans", .base = 2, .norm.entropy = F, .cor = "pearson", .do.norm = T, .laplace = 1e-12, .verbose = T, .k = 2, .eps = 0.01, .perp = 1, .theta = 0.1)
vis(gu.kmeans)
gu.kmeans = geneUsageAnalysis(tcr_gu, .method = "js+pca+kmeans", .base = 2, .norm.entropy = F, .cor = "spearman", .do.norm = T, .laplace = 1e-12, .verbose = T, .k = 2, .eps = 0.01, .perp = 1, .theta = 0.1)
vis(gu.kmeans)
s=vis(gu.kmeans)
fixVis(s)

gu.kmeans = geneUsageAnalysis(tcr_gu, .method = "js+pca+kmeans", .base = 2, .norm.entropy = F, .cor = "kendall", .do.norm = T, .laplace = 1e-12, .verbose = T, .k = 2, .eps = 0.01, .perp = 1, .theta = 0.1)
vis(gu.kmeans)

# Compare diversity of repertoires and visualise samples, grouped by two parameters
div = repDiversity(tcr$data, .method = "chao1")
vis(div, .by="pipeline", .meta=tcr$meta)
d=vis(div)
fixVis(d)
#vis(div, .by=c("Status", "Treatment"), .meta=tcr$meta)

entropy.seg(tcr$data, .genes=HUMAN_TRBV)
imm.js <-js.div.seg(tcr[1:10], HUMAN_TRBV, .verbose = F)
vis.radarlike(imm.js, .ncol = 2)

# Manipulate the visualisation of diversity estimates to make the plot publication-ready
#div.plot = vis(div, .by=c("Status", "Treatment"), .meta=immdata$meta)
#fixVis(div.plot)

trbv_gu = geneUsage(tcr$data, "hs.trbv", .norm = T, .ambig = "exc")
vis(trbv_gu, .plot = "hist", .grid = T)
trbj_gu = geneUsage(tcr$data, "hs.trbj", .norm = T, .ambig = "exc")
vis(trbj_gu, .plot = "hist", .grid = T)
trav_gu = geneUsage(tcr$data, "hs.trav", .norm = T, .ambig = "exc")
vis(trav_gu, .plot = "hist", .grid = T)
traj_gu = geneUsage(tcr$data, "hs.traj", .norm = T, .ambig = "exc")
vis(traj_gu, .plot = "hist", .grid = T)

trbv_gu = geneUsage(tcr$data, "hs.trbv", .norm = T, .ambig = "exc")
trbv_gu_cor = geneUsageAnalysis(trbv_gu, .method = "cor", .verbose = F)
vis(trbv_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=3.5, .signif.digits= 1)

v=vis(tcr_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=3.5, .signif.digits= 1)
fixVis(v)


trbj_gu = geneUsage(tcr$data, "hs.trbj", .norm = T, .ambig = "exc")
trbj_gu_cor = geneUsageAnalysis(trbj_gu, .method = "cor", .verbose = F)
v2 = vis(tcr_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=3.5, .signif.digits= 1)
fixVis(v2)

tcr_gu_js = geneUsageAnalysis(trbv_gu, .method = "js", .verbose = F)
vis(tcr_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size=3.5, .signif.digits=1)


gridExtra::grid.arrange(vis(tcr_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size=1.5), vis(tcr_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=1.5), ncol = 2)
