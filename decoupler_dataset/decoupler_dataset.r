library(Seurat)
library(decoupleR)
library(LR2TF)
inputs_dir <- system.file("extdata", package = "decoupleR")
data <- readRDS(file.path(inputs_dir, "sc_data.rds"))

#LR2TF::convert_seurat_to_anndata(data, "/home/larissa/Documents/LR2TF_HiWi/decoupler_test/")

p <- Seurat::DimPlot(data, 
                     reduction = "umap", 
                     label = TRUE, 
                     pt.size = 0.5) + 
     Seurat::NoLegend()

p

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

net

# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA@data)

# Run ulm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)

# Extract ulm and store it in tfsulm in pbmc
data[['tfsulm']] <- acts %>%
                    tidyr::pivot_wider(id_cols = 'source', 
                                       names_from = 'condition',
                                       values_from = 'score') %>%
                    tibble::column_to_rownames('source') %>%
                    Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfsulm"

# Scale the data
data <- Seurat::ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

p1 <- Seurat::DimPlot(data, 
                      reduction = "umap", 
                      label = TRUE, 
                      pt.size = 0.5) + 
      Seurat::NoLegend() + 
      ggplot2::ggtitle('Cell types')


colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p2 <- Seurat::FeaturePlot(data, features = c("PAX5")) + 
      ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
      ggplot2::ggtitle('PAX5 activity')


DefaultAssay(object = data) <- "RNA"
p3 <- Seurat::FeaturePlot(data, 
                          features = c("PAX5")) + 
      ggplot2::ggtitle('PAX5 expression')

Seurat::DefaultAssay(data) <- "tfsulm"

p <- p1 | p2 | p3
p

n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$tfsulm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster = Seurat::Idents(data)) %>%
      tidyr::pivot_longer(cols = -cluster, 
                          names_to = "source", 
                          values_to = "score") %>%
      dplyr::group_by(cluster, source) %>%
      dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(mean)) %>%
       dplyr::arrange(-abs(std)) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
                dplyr::filter(source %in% tfs) %>%
                tidyr::pivot_wider(id_cols = 'cluster', 
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames('cluster') %>%
                as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20) 