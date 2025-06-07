# R-scripts
Example scripts written in R for 
1. single-cell data analysis using Seurat package
2. analysis of TCR repertoire data using Immunarch package
3. survival analysis
4. RNA velocity esitmation using Velocyto.R package
5. Mutation mapper lollipop plots using cBioPortal and TCGA survival plots
6. Shiny App for single-cell data visualization

Single-Cell RNA-seq Explorer (Shiny App)

This interactive Shiny application allows users to visually explore and analyze single-cell RNA-seq data stored in a Seurat object. Key features include:

🧭 Main Functionalities:
Dimensionality Reduction Plots (DimPlot): Visualize cells using UMAP, t-SNE, or PCA, colored by various metadata (e.g., clusters, cell types).
Feature Plot: Plot gene expression across cells for selected genes.
Violin Plot: Visualize gene expression distributions across groups (e.g., clusters or cell types), with optional log transformation.
Heatmap: Generate heatmaps for selected gene sets, grouped by metadata, with options to customize display.
QC Tab: Includes:

Violin plots of quality control metrics (nFeature_RNA, nCount_RNA, percent.mt, percent.ribo)
Elbow Plot for PCA dimension selection
Mahalanobis Distance Plot to detect outlier cells based on QC metrics.
Metadata Table: View and search cell-level metadata in a sortable, scrollable table.

⚙️ Customizations & Enhancements:
Dynamically adjust point size, color grouping, and log-transformation.
Control gene selection via dropdown or comma-separated input.
Improved aesthetics: centered plots, consistent color palettes, larger fonts, and reduced label overlap.
