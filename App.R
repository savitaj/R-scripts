library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(patchwork)
library(shinythemes)
library(stats)
library(ComplexHeatmap)
library(tidyr)

setwd("D:/R-Scripts/New_ShinyApp")
seurat_data <- readRDS("seuratObj.rds")

# Make sure mitochondrial and ribosomal percent exist
if (!"percent.mt" %in% colnames(seurat_data@meta.data)) {
  seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^mt-")
}
if (!"percent.ribo" %in% colnames(seurat_data@meta.data)) {
  seurat_data[["percent.ribo"]] <- PercentageFeatureSet(seurat_data, pattern = "^Rps|^Rpl")
}

#DefaultAssay(seurat_data) <- "SCT"  # or "RNA", based on your pipeline
DefaultAssay(seurat_data) <- "RNA"

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Single-Cell RNA-seq Explorer"),
  sidebarPanel(
    h4("Plot Settings"),
    selectInput("reduction", "Dimensionality Reduction:",
                choices = c("umap", "tsne", "pca"),
                selected = "umap"),
    
    selectInput("group_by", "Group Cells By:",
                choices = c("ident", intersect(c("orig.ident", "samples", "seurat_clusters", "celltypes", "condition"),
                                               colnames(seurat_data@meta.data))),
                selected = "celltypes"),
    
    selectInput("feature", "Choose Gene Feature:",
                choices = rownames(seurat_data),
                selected = if ("Gapdh" %in% rownames(seurat_data)) "Gapdh" else rownames(seurat_data)[1]),
    
    sliderInput("pt.size", "Point Size:",
                min = 0.1, max = 3, value = 1, step = 0.1),
    
    # checkboxInput("log_transform", "Log Transform Feature Expression", TRUE),
    
    # tags$hr(),
    # 
    # h4("Heatmap Settings"),
    # textInput("heatmap_genes", "Genes (comma-separated):",
    #           value = "Cd3d,Cd8a,Ms4a1,Gapdh"),
    # 
    tags$hr(),
    
    h5("Tips:"),
    helpText("• Use gene symbols (case-sensitive)."),
    helpText("• Reduce point size for dense plots.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("QC Plots",
                 selectInput("qc_plot_type", "QC Plot Type:",
                             choices = c("Violin Plot", "Elbow Plot", "Mahalanobis Distance"),
                             selected = "Violin Plot"),
                 
                 conditionalPanel(
                   condition = "input.qc_plot_type == 'Violin Plot'",
                   helpText("Showing all QC metrics (nFeature, nCount, percent.mt, percent.ribo) per sample.")
                 ),
                 
                 plotOutput("qc_plot", height = "600px")
        ),
        tabPanel("Dim Plot", plotOutput("dimplot")),
        tabPanel("Feature Plot", plotOutput("featureplot")),
        tabPanel("Violin Plot", plotOutput("vlnplot")),
        tabPanel("Heatmap",
                 fluidRow(
                   column(
                     width = 12,
                     textInput("heatmap_genes", "Heatmap Genes (comma-separated):", value = "Cd3d,Cd8a,Ms4a1,Gapdh")
                   )
                 ),
                 fluidRow(
                   column(
                     width = 12,
                     align = "center",
                     plotOutput("heatmap")
                   )
                 )
        ),
        #tabPanel("Heatmap", plotOutput("heatmap")),
        tabPanel("Meta Data", DTOutput("metadata_table"))
      )
    )
)

server <- function(input, output, session) {

  output$dimplot <- renderPlot({
    DimPlot(seurat_data, reduction = input$reduction, group.by = input$group_by, pt.size = input$pt.size)
  })
  
  output$featureplot <- renderPlot({
    FeaturePlot(
      object = seurat_data,
      features = input$feature,
      reduction = input$reduction,
      pt.size = input$pt.size
    )
  })
  
  output$vlnplot <- renderPlot({
    # Fetch expression for the selected gene
    data_expr <- FetchData(seurat_data, vars = input$feature)
    
    # Apply log1p transform if selected
    # if (input$log_transform) {
    #   data_expr[[input$feature]] <- log1p(data_expr[[input$feature]])
    # }
    
    # Grouping variable from metadata
    group_var <- seurat_data@meta.data[[input$group_by]]
    data_expr$group <- factor(group_var)
    
    # Plot
    ggplot(data_expr, aes(x = group, y = .data[[input$feature]], fill = group)) +
      geom_violin(trim = FALSE) +
      theme_minimal(base_size = 16) +
      xlab(input$group_by) + ylab("Expression") +
      ggtitle(paste("Violin Plot of", input$feature)) +
      scale_fill_manual(values = scales::hue_pal()(length(levels(data_expr$group)))) +
      theme(
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold")
      )
  })
  
  # output$heatmap <- renderPlot({
  #   gene_list <- unlist(strsplit(input$heatmap_genes, ",\\s*"))
  #   gene_list <- gene_list[gene_list %in% rownames(seurat_data)]
  #   if (length(gene_list) > 1) {
  #     expr <- seurat_data[["RNA"]]@scale.data[gene_list, ]
  #     Heatmap(expr,
  #             show_column_names = FALSE,
  #             cluster_columns = FALSE,
  #             name = "Expression")
  #   } else {
  #     plot.new()
  #     text(0.5, 0.5, "Please enter 2 or more valid genes.")
  #   }
  # })


  
  output$heatmap <- renderPlot({
    gene_list <- unlist(strsplit(input$heatmap_genes, ",\\s*"))
    gene_list <- gene_list[gene_list %in% rownames(seurat_data)]
    
    if (length(gene_list) > 1) {
      p <- DoHeatmap(
        seurat_data,
        features = gene_list,
        group.by = input$group_by
      ) +
        theme(
          axis.text.x = element_blank(),             # Remove top labels
          axis.text.y = element_text(size = 8),      # Smaller gene labels
          plot.margin = margin(20, 10, 10, 10),      # Move plot downward
          legend.position = "bottom"                 # Move legend below
        ) +
        guides(fill = guide_colorbar(title.position = "top"))
      
      p
    } else {
      plot.new()
      text(0.5, 0.5, "Please enter 2 or more valid genes.")
    }
  })
  
  
  
  
  output$qc_plot <- renderPlot({
    req(seurat_data)  # Ensure object is loaded
    
    if (input$qc_plot_type == "Violin Plot") {
      qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
      
      # Check if each metric exists in meta.data
      available <- qc_metrics[qc_metrics %in% colnames(seurat_data@meta.data)]
      
      if (length(available) == 0) {
        plot.new()
        text(0.5, 0.5, "No QC metrics available in meta.data.")
        return()
      }
      
      df <- seurat_data@meta.data
      df$sample <- df$orig.ident
      
      # Ensure the QC columns aren't all NA
      valid_metrics <- available[sapply(df[available], function(col) !all(is.na(col)))]
      
      if (length(valid_metrics) == 0) {
        plot.new()
        text(0.5, 0.5, "All QC metrics are empty or NA.")
        return()
      }
      
      df_long <- tidyr::pivot_longer(df,
                                     cols = all_of(valid_metrics),
                                     names_to = "metric",
                                     values_to = "value")
      
      ggplot(df_long, aes(x = sample, y = value, fill = sample)) +
        geom_violin(trim = FALSE) +
        facet_wrap(~ metric, scales = "free_y", ncol = 2) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(face = "bold")) +
        xlab("Sample") + ylab("Value") +
        ggtitle("QC Metrics per Sample")
      
    } else if (input$qc_plot_type == "Elbow Plot") {
      ElbowPlot(seurat_data, ndims = 50)
      
    } else if (plot_type == "Mahalanobis Distance") {
      qc_df <- seurat_data@meta.data %>%
        dplyr::select(nFeature_RNA, nCount_RNA, percent.mt, percent.ribo) %>%
        na.omit()
      
      mdist <- mahalanobis(qc_df, colMeans(qc_df), cov(qc_df))
      
      ggplot(data.frame(MahalanobisDistance = mdist),
             aes(x = MahalanobisDistance)) +
        geom_histogram(fill = "tomato", bins = 50, color = "white") +
        theme_minimal() +
        ggtitle("Mahalanobis Distance QC") +
        xlab("Mahalanobis Distance")
    }
  })
  
  output$metadata_table <- renderDT({
    datatable(seurat_data@meta.data, options = list(scrollX = TRUE))
  })
}

shinyApp(ui = ui, server = server)
