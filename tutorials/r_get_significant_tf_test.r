
library(LR2TF)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(maditr)

TF_object <- readRDS("LR2TF_test_run/results/TF_results/result_TF_object.RDS")
TF_object@intracellular_network_condition
TF_object@intracellular_network_cluster



data(bone_marrow_stromal_cell_example, package = "LR2TF")
seuratobject <- bone_marrow_stromal_cell_example

seuratobject@meta.data

DefaultAssay(object = seuratobject) <- "tf_activities"
seuratobject <- ScaleData(seuratobject)

Idents(object = seuratobject) <- protocol
seuratobject[['tf_condition']] <- Idents(object = seuratobject)
Idents(object = seuratobject) <- new_annotation
seuratobject[['tf_annotation']] <- Idents(object = seuratobject)

levels(seuratobject@meta.data$tf_annotation)

pval <- 0.05
log2fc <- 0
tf_activities <- t(read.csv("decoupler_results.csv", header = TRUE, row.names = 1))
tf_acticities <- CreateAssayObject(data = tf_activities)
seuratobject[["tf_activities"]] <- tf_acticities

DefaultAssay(object = seuratobject) <- "tf_activities"
seuratobject <- ScaleData(seuratobject)
seuratobject[['tf_annotation']] <- Idents(object = seuratobject)
I

Idents(object = seuratobject) <- "protocol"
seuratobject_list <- SplitObject(seuratobject, split.by ="protocol")
seuratobject <- seuratobject_list[[1]]

Idents(object = seuratobject) <- "tf_annotation"

  number_of_clusters = length(levels(Idents(seuratobject)))

  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE,
                                         min.pct = 0, logfc.threshold = 0,
                                         verbose = FALSE)

  seuratobject.markers$tag <- sapply(seuratobject.markers$p_val_adj, function(pval) {
    if (pval < 0.001) {
      txt <- "***"
    } else if (pval < 0.01) {
      txt <- "**"
    } else if (pval < 0.05) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  seuratobject.markers$log_fc_tag <- sapply(seuratobject.markers$avg_log2FC, function(log_fc) {
    if (log_fc >= 1.0) {
      txt <- "***"
    } else if (log_fc > 0.5) {
      txt <- "**"
    } else if (log_fc > 0.0) {
      txt <- "*"
    } else {
      txt <- "ns"
    }
    return(txt)
  })

  tag_mapping = seuratobject.markers[c("gene", "tag", "log_fc_tag", "cluster", "avg_log2FC", "p_val_adj")]
  tag_mapping = filter(tag_mapping, p_val_adj < as.double(pval))
  tag_mapping = filter(tag_mapping, avg_log2FC > as.double(log2fc) | avg_log2FC < (0 - as.double(log2fc)))
  #tag_mapping = tag_mapping[(tag_mapping$tag == "***"),]
  #tag_mapping = tag_mapping[(tag_mapping$log_fc_tag == "***"),]
  tag_mapping = dcast(tag_mapping, gene ~ cluster, value.var = "tag")
  row.names(tag_mapping) = tag_mapping$gene
  tag_mapping$gene = NULL
  tag_mapping$log_fc_tag = NULL
  tag_mapping[is.na(tag_mapping)] <- "ns"

  viper_scores_df <- GetAssayData(seuratobject, slot = "scale.data",
                                  assay = "tf_activities") %>%
    data.frame(check.names = F) %>%
    t()

  CellsClusters <- data.frame(cell = names(Idents(seuratobject)),
                              cell_type = as.character(Idents(seuratobject)),
                              check.names = F)



  #col.num <- which(colnames(viper_scores_df) %in% rownames(tag_mapping))
  #viper_scores_df <- viper_scores_df[, sort(c(col.num))]

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
print(summarized_viper_scores, n=300)

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df = t(summarized_viper_scores_df)
  write.csv(summarized_viper_scores_df, file = paste0('testtttt.csv'))


  highly_variable_tfs_all <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg)) %>%
    ungroup() %>%
    distinct(tf)

  summarized_viper_scores_df_all <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs_all, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)
  tf_scores <- t(summarized_viper_scores_df_all)
  write.csv(tf_scores, file = paste0('variable_tf_test.csv'))
  
test <- c(-0.116, 1.06, -0.783, -0.759, -0.0925)
test_var <- var(test)
#even this is different from the var(avg) in the variable tf function


############################################################################################################
ligands <- load(file='ligands_human.rda')
str(ligands)
