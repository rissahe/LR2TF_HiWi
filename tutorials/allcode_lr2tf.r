#' Generate cluster and condition heatmap with r effect size only for significant genes
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param celltype_annotation meta data field with celltype annotations
#' @param condition_annotation meta data field with condition annotation
#' @param comparison_list list of wished comparisons
#' @import glue
#' @import maditr
#' @export
condition_comparison_significant <- function(seuratobject, out_path, celltype_annotation, condition_annotation, comparison_list) {

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- condition_annotation
  seuratobject[['tf_condition']] <- Idents(object = seuratobject)
  Idents(object = seuratobject) <- celltype_annotation
  seuratobject[['tf_annotation']] <- Idents(object = seuratobject)

  vs_df_list <- list()

  for (vs in comparison_list) {
    vs1 <- vs[1]
    vs2 <- vs[2]

    message("vs: ", vs1, " ", vs2, " ", date(), "\n")
    Idents(seuratobject) <- condition_annotation

    pws <- rownames(seuratobject@assays$tf_activities)

    ###
    res <- list()
    for (i in levels(seuratobject@meta.data$tf_annotation)) {
      a_sub <- subset(seuratobject, cells = rownames(seuratobject@meta.data)[seuratobject@meta.data$tf_annotation == i & (seuratobject@meta.data$tf_condition %in% vs)])
      if (length(unique(a_sub@meta.data$tf_condition)) == 2) {
        condition_table <- as.data.frame(a_sub@meta.data$tf_condition)
        names(condition_table)[1] <- "condition"
        metadata_counts <- condition_table %>%
          group_by(condition) %>%
          summarise(total_count = n())
        if (all(metadata_counts$total_count > 10)) {
          g <- as.character(a_sub@meta.data$tf_condition)
          g <- factor(g, levels = c(vs1, vs2)) ###############################################
          res[[i]] <- scran::findMarkers(as.matrix(a_sub@assays$tf_activities@scale.data), g)[[1]]
          res[[i]] <- as.data.frame(res[[i]])
          r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$tf_activities@scale.data[pw,]), g))
          nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
          names(r) <- nms
          res[[i]][nms, "r"] <- r
          res[[i]] <- res[[i]][nms,]
        }
      }
    }

    for (cl in names(res)) {
      res[[cl]]$tf <- rownames(res[[cl]])
      res[[cl]]$CellType <- cl
      colnames(res[[cl]]) <- c("Top", "p.value", "FDR", "summary.logFC", "logFC", "r", "tf", "CellType")

    }
    res_df <- do.call("rbind", res)
    res_df <- na.omit(res_df)
    res_df$tag <- sapply(res_df$FDR, function(pval) {
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

    #significant_res <- res_df[res_df$tag == '***',]
    #significant_genes <- unique(significant_res$tf)

    #end_res <- filter(res_df, res_df$tf %in% significant_genes)
    end_res <- res_df

    vs_df_list[[glue("{vs1} vs {vs2}")]] <- end_res
    write.csv(res_df, paste0(out_path, "/all_tfs_", glue("{vs1}_vs_{vs2}", ".csv")))
  }

  saveRDS(vs_df_list, file = paste0(out_path, "/comparison_dfs.RDS"))
  return(vs_df_list)
}

#' Mouse-Human homologous genes table transcription factors
#'
#' A data frame containing the mouse-human homologous genes of the dorothea regulon
#'
#' @format A data frame with 1026 rows and 4 variables:
#' \itemize{
#'   \item 10090: Mouse gene symbol
#'   \item 9606: Human gene symbol
#'   \item 10090_ID: Mouse NCBI gene ID
#'   \item 9606_ID: Human NCBI gene ID
#' }
"mouse_regulon_homologous"

#' List of annotated ligands from OmnipathR for human data
#'
#' @format A list with 9708 elements
"ligands_human"

#' Mouse-Human homologous genes table ligands
#'
#' A data frame containing the mouse-human homologous genes of the ligands in the Omnipath database
#'
#' @format A data frame with 1026 rows and 4 variables:
#' \itemize{
#'   \item 10090: Mouse gene symbol
#'   \item 9606: Human gene symbol
#' }
"ligands_homologous"

#' Example dataset
#'
#' Example dataset for tutorial execution of the package
#'
#'
"bone_marrow_stromal_cell_example"

#' Dataframe with pre-computed receptor to transcription factor connections
#'
#' Dataframe with pre-computed receptor to transcription factor connections from Omnipath
#'
"RTF_DB"

#' Dataframe with pre-computed receptor to transcription factor connections 2nd version
#'
#' Dataframe with pre-computed receptor to transcription factor connections from Omnipath
#'
"RTF_DB_2"

#' Run DoRothEA analysis to get significant transcription factors
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param condition Experminet condition (e.g. disease, knockout ...)
#' @param out_path Output path to save results
#' @param tf_condition_significant condition comparison results
#' @param pval p-value to filter results
#' @param log2fc log fold change value to filter results
#' @return A data frame with transcription factor activity scores per cell type
#' @export
get_significant_tfs <- function(seuratobject, condition, out_path, tf_condition_significant, pval, log2fc) {
  single_result_path <- paste0(out_path, condition)
  dir.create(single_result_path)

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- "tf_annotation"

  number_of_clusters = length(levels(Idents(seuratobject)))

  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE,
                                         min.pct = 0, logfc.threshold = 0,
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file =
    paste0(single_result_path, '/all_specificmarker_',
           '_', condition, '.csv'))

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

  save_unfiltered_tf_scores(viper_scores_df, CellsClusters, condition, single_result_path)

  unfiltered_viper_scores = unfiltered_tf_activity_table(viper_scores_df, CellsClusters)

  col.num <- which(colnames(viper_scores_df) %in% rownames(tag_mapping))
  viper_scores_df <- viper_scores_df[, sort(c(col.num))]

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df = t(summarized_viper_scores_df)
  write.csv(summarized_viper_scores_df, file = paste0(single_result_path, '/tf_scores', '_', condition, '.csv'))

  tf_scores_variable_table = save_variable_tf_scores(summarized_viper_scores, condition, single_result_path)

  message("Plotting top variables tf activities")
  plot_highly_variable_tfs(summarized_viper_scores, condition,
                           single_result_path, number_of_clusters)

  plot_tf_activity_compressed(summarized_viper_scores_df, condition, single_result_path)
  plot_tf_activity(summarized_viper_scores_df, tag_mapping, condition,
                   single_result_path)

  rownames(summarized_viper_scores_df) = gsub(".", "-", rownames(summarized_viper_scores_df), fixed = TRUE)
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$tag == "***"),]
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$log_fc_tag == "***"),]

  map_z_value_filtered <- function(gene, cluster) {
    if (gene %in% rownames(summarized_viper_scores_df)) {
      z_score = summarized_viper_scores_df[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  map_z_value <- function(gene, cluster) {
    if (gene %in% rownames(unfiltered_viper_scores)) {
      z_score = unfiltered_viper_scores[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  res <- list()
  seuratobject.markers$z_score = mapply(map_z_value_filtered, seuratobject.markers$gene, seuratobject.markers$cluster)
  seuratobject.markers = seuratobject.markers[!(seuratobject.markers$tag == "ns"),]
  seuratobject.markers = na.omit(seuratobject.markers)
  res[["cluster"]] = seuratobject.markers[c("gene", "tag", "cluster", "z_score")]
  write.csv(res[["cluster"]], file = paste0(single_result_path, '/significant_cluster_tf_results', '_', condition, '.csv'))

  tf_condition_significant$gene = gsub(".", "-", tf_condition_significant$gene, fixed = TRUE)
  tf_condition_significant$z_score = mapply(map_z_value, tf_condition_significant$gene, tf_condition_significant$cluster)
  tf_condition_significant = na.omit(tf_condition_significant)
  res[["condition"]] = tf_condition_significant
  write.csv(res[["condition"]], file = paste0(single_result_path, '/significant_condition_tf_results', '_', condition, '.csv'))

  return(res)
}

#' Run DoRothEA analysis to get significant transcription factors
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param condition Experminet condition (e.g. disease, knockout ...)
#' @param out_path Output path to save results
#' @param tf_condition_significant condition comparison results
#' @return A data frame with transcription factor activity scores per cell type
#' @export
get_significant_tfs_single <- function(seuratobject, condition, out_path, pval, log2fc) {
  single_result_path <- paste0(out_path, condition)
  dir.create(single_result_path)

  DefaultAssay(object = seuratobject) <- "tf_activities"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- "tf_annotation"

  number_of_clusters = length(levels(Idents(seuratobject)))

  seuratobject.markers <- FindAllMarkers(seuratobject, only.pos = TRUE,
                                         min.pct = 0, logfc.threshold = 0,
                                         verbose = FALSE)
  write.csv(seuratobject.markers, file =
    paste0(single_result_path, '/all_specificmarker_',
           '_', condition, '.csv'))

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

  tag_mapping = seuratobject.markers[c("gene", "tag", "log_fc_tag", "cluster", "p_val_adj", "avg_log2FC")]
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

  save_unfiltered_tf_scores(viper_scores_df, CellsClusters, condition, single_result_path)

  col.num <- which(colnames(viper_scores_df) %in% rownames(tag_mapping))
  viper_scores_df <- viper_scores_df[, sort(c(col.num))]

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df = t(summarized_viper_scores_df)
  write.csv(summarized_viper_scores_df, file = paste0(single_result_path, '/tf_scores', '_', condition, '.csv'))

  tf_scores_variable_table = save_variable_tf_scores(summarized_viper_scores, condition, single_result_path)

  message("Plotting top variables tf activities")
  plot_highly_variable_tfs(summarized_viper_scores, condition,
                           single_result_path, number_of_clusters)

  plot_tf_activity_compressed(summarized_viper_scores_df, condition, single_result_path)
  plot_tf_activity(summarized_viper_scores_df, tag_mapping, condition,
                   single_result_path)

  rownames(summarized_viper_scores_df) = gsub(".", "-", rownames(summarized_viper_scores_df), fixed = TRUE)
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$tag == "***"),]
  #seuratobject.markers = seuratobject.markers[(seuratobject.markers$log_fc_tag == "***"),]

  map_z_value_filtered <- function(gene, cluster) {
    if (gene %in% rownames(summarized_viper_scores_df)) {
      z_score = summarized_viper_scores_df[as.character(gene), as.character(cluster)]
      return(z_score)
    } else {
      return(NA)
    }
  }

  res <- list()
  seuratobject.markers$z_score = mapply(map_z_value_filtered, seuratobject.markers$gene, seuratobject.markers$cluster)
  seuratobject.markers = seuratobject.markers[!(seuratobject.markers$tag == "ns"),]
  seuratobject.markers = na.omit(seuratobject.markers)
  res[["cluster"]] = seuratobject.markers[c("gene", "tag", "cluster", "z_score")]
  write.csv(res[["cluster"]], file = paste0(single_result_path, '/significant_cluster_tf_results', '_', condition, '.csv'))

  return(res)
}

#' Run Dorothea based on Seurat object
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param confidence_level Curation confidence level to filter DoRothEA regulon (default "ABC")
#' @param organism Organism of sample origin
#' @return Seurat object with dorothea assay
#' @export
dorothea_base_execution <- function(seuratobject, out_path, confidence_level = c("A", "B", "C"),
                                    organism = "human") {

  if (organism == "human") {
    dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% confidence_level)
  } else if (organism == "mouse") {
    dorothea_regulon_human <- get(data("dorothea_mm", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% confidence_level)
  } else {
    dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% confidence_level)
  }

  seuratobject <- run_viper(seuratobject, regulon,
                            options = list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))

  saveRDS(seuratobject, file = paste0(out_path, "seuratobject_dorothea_results.RDS"))

  return(seuratobject)
}


#' Generate cluster and condition heatmap with r effect size only for significant genes
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @param celltype_annotation meta data field with celltype annotations
#' @param condition_annotation meta data field with condition annotation
#' @param comparison_list list of wished comparisons
#' @param organism Organism of sample origin
#' @return tables with r-effectsize results
#' @export
calculate_r_effectsize <- function(seuratobject, organism, out_path, celltype_annotation, condition_annotation, comparison_list) {

  DefaultAssay(object = seuratobject) <- "dorothea"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- condition_annotation
  seuratobject[['doro_condition']] <- Idents(object = seuratobject)
  Idents(object = seuratobject) <- celltype_annotation
  seuratobject[['tf_annotation']] <- Idents(object = seuratobject)

  vs_df_list <- list()

  for (vs in comparison_list) {
    vs1 <- vs[1]
    vs2 <- vs[2]

    message("vs: ", vs1, " ", vs2, " ", date(), "\n")
    Idents(seuratobject) <- condition_annotation

    pws <- rownames(seuratobject@assays$dorothea)

    ###
    res <- list()
    for (i in levels(seuratobject@meta.data$tf_annotation)) {
      a_sub <- subset(seuratobject, cells = rownames(seuratobject@meta.data)[seuratobject@meta.data$tf_annotation == i & (seuratobject@meta.data$doro_condition %in% vs)])
      g <- as.character(a_sub@meta.data$doro_condition)
      g <- factor(g, levels = c(vs1, vs2)) ###############################################
      res[[i]] <- scran::findMarkers(as.matrix(a_sub@assays$dorothea@scale.data), g)[[1]]
      res[[i]] <- as.data.frame(res[[i]])
      r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$dorothea@scale.data[pw,]), g))
      nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
      names(r) <- nms
      res[[i]][nms, "r"] <- r
      res[[i]] <- res[[i]][nms,]

    }

    for (cl in names(res)) {
      res[[cl]]$tf <- rownames(res[[cl]])
      res[[cl]]$CellType <- cl
      colnames(res[[cl]]) <- c("Top", "p.value", "FDR", "summary.logFC", "logFC", "r", "tf", "CellType")

    }
    res_df <- do.call("rbind", res)
    res_df$tag <- sapply(res_df$FDR, function(pval) {
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

    vs_df_list[[glue("{vs1} vs {vs2}")]] <- res_df
    write.csv(res_df, paste0(out_path, "all_tfs_", glue("{vs1}_vs_{vs2}", ".csv")))
  }
  return(vs_df_list)
}


#' Group transcription factor activity by gene and cluster.
#'
#' Description
#'
#' @param activity_data_frame Data frame with tf activity score by cell
#' @param CellsClusters Cluster annotation by cell
#' @return A data frame with transcription factor activity scores per cell type
#' @export
unfiltered_tf_activity_table <- function(activity_data_frame, CellsClusters) {
  viper_scores_clusters <- activity_data_frame %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  summarized_viper_scores_df <- t(summarized_viper_scores_df)

  return(summarized_viper_scores_df)
}


#' Saves transcription factor activity scores per cell type into a table
#'
#' This function saves the transcription factor activity scores per cell type
#' into a csv table.
#'
#' @param viper_scores_df data frame with transcription factor activity scores per cell
#' @param CellsClusters Clusters perr cell
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import stringr
#' @export
save_unfiltered_tf_scores <- function(viper_scores_df, CellsClusters, condition, out_path) {

  viper_scores_clusters <- viper_scores_df %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_df <- summarized_viper_scores %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)
  tf_scores <- t(summarized_viper_scores_df)
  write.csv(tf_scores, file = paste0(out_path, '/unfiltered_tf_scores', '_', condition, '.csv'))
}


#' Saves most variable transcription factor activity scores per cell type into a table
#'
#' This function saves the transcription factor activity scores per cell type
#' into a csv table.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import stringr
#' @export
save_variable_tf_scores <- function(tf_scores, condition, out_path) {

  highly_variable_tfs_all <- tf_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg)) %>%
    ungroup() %>%
    distinct(tf)

  summarized_viper_scores_df_all <- tf_scores %>%
    semi_join(highly_variable_tfs_all, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)
  tf_scores <- t(summarized_viper_scores_df_all)
  write.csv(tf_scores, file = paste0(out_path, '/variable_tf_scores', '_', condition, '.csv'))

  return(tf_scores)
}


#' Plot heatmap of all transcription factor activities
#'
#' This function plots a heatmap with the transcription factor activities per
#' cell type.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param tag_mapping Labels for heatmap
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import pheatmap
#' @import tidyr
#' @import stringr
#' @export
plot_tf_activity <-
  function(tf_scores,
           tag_mapping,
           condition,
           out_path) {
    tf_scores <- as.matrix(tf_scores)

    plot_width <- ((ncol(tf_scores) * 15) / 25.4) + 5
    plot_height <- ((nrow(tf_scores) * 3) / 25.4) + 5

    column_names_table <- colnames(tf_scores)
    column_names_tags <- colnames(tag_mapping)
    if (length(column_names_tags) < length(column_names_table)) {
      missing_col <- setdiff(column_names_table, column_names_tags)
      for (col in missing_col) {
        tag_mapping[col] <- "ns"
      }
    }

    rownames(tf_scores) <- gsub(".", "-", rownames(tf_scores), fixed = TRUE)

    pdf(
      file = paste0(out_path, '/tf_activity_', condition, '.pdf'),
      height = plot_height,
      width = plot_width
    )

    fh <- function(x)
      fastcluster::hclust(dist(x))

    p <-
      Heatmap(
        tf_scores,
        name = "z-score",
        cluster_rows = fh,
        width = ncol(tf_scores) * unit(15, "mm"),
        height = nrow(tf_scores) * unit(3, "mm"),
        row_title = "Transcription Factor",
        column_title = "Cell Type",
        row_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(tag_mapping[as.character(rownames(tf_scores)[i]), as.character(colnames(tf_scores)[j])], x, y, gp = gpar(fontsize = 8))
        }
        # layer_fun = function(j, i, x, y, width, height, fill) {
        #   v = pindex(tf_scores, i, j)
        #   grid.text(tag_mapping[as.character(rownames(tf_scores)[i]),as.character(colnames(tf_scores)[j])], x, y, gp = gpar(fontsize = 10))
        # }
      )
    draw(p)
    dev.off()
  }


#' Adds gene type to the name of the gene
#'
#' Description
#'
#' @param df dataframe with all interactions
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
add_node_type <- function(df) {
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Ligand", paste0(gene_A, "|L"), gene_A))
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Receptor", paste0(gene_A, "|R"), gene_A))
  df <- df %>%
    mutate(gene_A = ifelse(type_gene_A == "Transcription Factor", paste0(gene_A, "|TF"), gene_A))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Ligand", paste0(gene_B, "|L"), gene_B))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Receptor", paste0(gene_B, "|R"), gene_B))
  df <- df %>%
    mutate(gene_B = ifelse(type_gene_B == "Transcription Factor", paste0(gene_B, "|TF"), gene_B))
}


#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_prediction path to or dataframe with ligand-receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
combine_LR_and_TF <- function(tf_table, LR_prediction, out_path, condition, add_node_type = FALSE) {

  if (!is.data.frame(LR_prediction)) {
    lr_table <- read.csv(LR_prediction)
    row.names(lr_table) <- lr_table$X
    lr_table$X <- NULL
  } else {
    lr_table <- LR_prediction
  }

  lr_ligands <- unique(lr_table$gene_A)
  lr_receptors <- unique(lr_table$gene_B)
  tf_receptor_interactions <- tf_table %>%
    filter(gene_A %in% lr_receptors)
  tf_ligand_interactions <- tf_table %>%
    filter(gene_B %in% lr_ligands)

  complete_interactions <- rbind(tf_receptor_interactions, tf_ligand_interactions)
  complete_interactions <- rbind(complete_interactions, lr_table)
  if (add_node_type) {
    complete_interactions <- add_node_type(complete_interactions)
  }
  write.csv(complete_interactions, paste0(out_path, "CrossTalkeR_input_", condition, ".csv"), row.names = FALSE)
  return(complete_interactions)
}

#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_path path to ligand receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
combine_LR_and_TF_unfiltered <- function(tf_table, LR_path, out_path, condition) {

  if (!is.data.frame(LR_prediction)) {
    lr_table <- read.csv(LR_prediction)
    row.names(lr_table) <- lr_table$X
    lr_table$X <- NULL
  } else {
    lr_table <- LR_prediction
  }

  complete_interactions <- rbind(tf_table, lr_table)
  complete_interactions <- add_node_type(complete_interactions)

  write.csv(complete_interactions, paste0(out_path, "CrossTalkeR_input_", condition, ".csv"), row.names = FALSE)
}


#' Extract all condition specific tfs in one table
#'
#' Description
#'
#' @param condition_tfs list of tables for condition specific tfs
#' @param out_path save path for result table
#' @param log2fc log2 fold-change threshold
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
unique_condition_specific_tfs <- function(condition_tfs, out_path, log2fc) {

  compared_tfs <- data.frame(
    gene = character(),
    tag = character(),
    cluster = character()
  )

  for (result_name in names(condition_tfs)) {

    tf_condition_significant <- condition_tfs[[result_name]]
    tf_condition_significant <- tf_condition_significant %>%
      filter(tag == "***") %>%
      filter(logFC > as.double(log2fc) | logFC < (0 - as.double(log2fc))) %>%
      subset(select = c("tf", "tag", "CellType")) %>%
      rename(gene = tf, cluster = CellType)

    compared_tfs <- rbind(compared_tfs, tf_condition_significant)
  }

  compared_tfs <- compared_tfs[!duplicated(compared_tfs$gene),]
  write.csv(compared_tfs, paste0(out_path, "all_condition_specific_tfs.csv"), row.names = FALSE)

  return(compared_tfs)
}


#' Plotting heatmaps for tf activity comparison between conditions
#'
#' Description
#'
#' @param seuratobject seurat object with drothea assay
#' @param condition_annotation seurat object meta data field name of conditions
#' @param celltype_annotation seurat object meta data field name of cluster
#' @param filter_list list of tf genes to filter all detected tf activities
#' @param out_path output path for the results
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
plot_condition_Heatmaps <- function(seuratobject, condition_annotation, celltype_annotation, filter_list, out_path) {

  tmp_out_path <- paste0(out_path, "tmp")
  dir.create(tmp_out_path)

  DefaultAssay(object = seuratobject) <- "dorothea"
  seuratobject <- ScaleData(seuratobject)

  Idents(object = seuratobject) <- condition_annotation
  seuratobject[['doro_condition']] <- Idents(object = seuratobject)
  Idents(object = seuratobject) <- celltype_annotation
  seuratobject[['tf_annotation']] <- Idents(object = seuratobject)

  seuratobject_list <- SplitObject(seuratobject, split.by = "tf_annotation")

  for (name in names(seuratobject_list)) {

    sub_object <- seuratobject_list[[name]]
    Idents(object = sub_object) <- "doro_condition"

    number_of_conditions <- length(levels(Idents(sub_object)))

    viper_scores_df <- GetAssayData(sub_object, slot = "scale.data",
                                    assay = "dorothea") %>%
      data.frame(check.names = F) %>%
      t()

    CellsClusters <- data.frame(cell = names(Idents(sub_object)),
                                cell_type = as.character(Idents(sub_object)),
                                check.names = F)

    viper_scores_clusters <- viper_scores_df %>%
      data.frame() %>%
      rownames_to_column("cell") %>%
      gather(tf, activity, -cell) %>%
      inner_join(CellsClusters)

    summarized_viper_scores <- viper_scores_clusters %>%
      group_by(tf, cell_type) %>%
      summarise(avg = mean(activity),
                std = sd(activity))

    plot_highly_variable_tfs_condition(summarized_viper_scores, name, tmp_out_path, number_of_conditions)

    summarized_viper_scores_df <- summarized_viper_scores %>%
      filter(tf %in% filter_list) %>%
      dplyr::select(-std) %>%
      spread(tf, avg) %>%
      data.frame(row.names = 1, check.names = FALSE)

    summarized_viper_scores_df <- t(summarized_viper_scores_df)

    plot_tf_activity_without_mapping_condition(summarized_viper_scores_df, name, tmp_out_path)
  }
  file_list <- list.files(path = tmp_out_path, pattern = "*.pdf", full.names = TRUE)
  qpdf::pdf_combine(input = file_list,
                    output = paste0(out_path, "condition_tf_heatmaps.pdf"))

  unlink(tmp_out_path, recursive = TRUE)
}

#' Create an empty dataframe with CrossTalkeR input table format
#'
#' @export
create_empty_CTR_dataframe <- function() {
  empty_df <- data.frame(
    source = character(),
    target = character(),
    gene_A = character(),
    gene_B = character(),
    type_gene_A = character(),
    type_gene_B = character(),
    MeanLR = numeric()
  )
  return(empty_df)
}

#' Add entry to a dataframe in the CrossTalkeR input format
#'
#' @export
add_entry_to_CTR_dataframe <- function(source, target, gene_A, gene_B, type_gene_A, type_gene_B, MeanLR) {
  df <-
    data.frame(
      source,
      target,
      gene_A,
      gene_B,
      type_gene_A,
      type_gene_B,
      MeanLR
    )
  names(df) <-
    c(
      'source',
      'target',
      'gene_A',
      'gene_B',
      'type_gene_A',
      'type_gene_B',
      "MeanLR"
    )

  return(df)
}


#' Create an empty dataframe with CrossTalkeR input table format
#'
#' @export
create_empty_Regulon_dataframe <- function() {
  empty_df <- data.frame(
    celltype = character(),
    Receptor = character(),
    TF = character(),
    Target_Gene = character(),
    TF_Score = numeric()
  )
  return(empty_df)
}

#' Add entry to a dataframe in the CrossTalkeR input format
#'
#' @export
add_entry_to_Regulon_dataframe <- function(celltype, Receptor, TF, Target_Gene, TF_Score) {
  df <-
    data.frame(
      celltype,
      Receptor,
      TF,
      Target_Gene,
      TF_Score
    )
  names(df) <-
    c(
      'celltype',
      'Receptor',
      'TF',
      'Target_Gene',
      'TF_Score'
    )

  return(df)
}

#' Convert Seurat object to anndata object and save anndata object file
#'
#' @param seuratobject Input Seurat Object
#' @param out_path Output path to save results
#' @import sceasy
#' @export
convert_seurat_to_anndata <- function(seuratobject, out_path) {
  sceasy::convertFormat(seuratobject, from = "seurat", to = "anndata", outFile = paste0(out_path, "anndata_object.h5ad"))
}


#' Run viper analysis with the decoupleR package
#'
#' @param seuratobject Input Seurat Object
#' @param regulon table with weighted TF to target gene interactions
#' @return seuratobject
#' #' @import decoupleR
#' @export
decoupleR_viper_analysis <- function(seuratobject, regulon) {
  mat <- as.matrix(seuratobject@assays$RNA@data)
  tf_activities <- decoupleR::run_viper(mat = mat, network = regulon, .source = 'source', .target = 'target', .mor = 'weight', verbose = TRUE)
  seuratobject[['tf_activities']] <- tf_activities %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

  return(seuratobject)
}


#' Loading old dorothea package regulons (default if no regulon is provided)
#'
#' @param organism of input data
#' @return regulon
#' @import dorothea
#' @export
load_dorothea_regulon <- function(organism) {
  if (organism == "human") {
    dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>%
      dplyr::filter(confidence %in% c("A", "B", "C", "D")) %>%
      rename(source = tf, weight = mor)
  }
  else if (organism == "mouse") {
    dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
    regulon <- dorothea_regulon_mouse %>%
      dplyr::filter(confidence %in% c("A", "B", "C", "D")) %>%
      rename(source = tf, weight = mor)
  } else {
    print("only human and mouse regulons can be loaded by default!")
  }
  return(regulon)
}

#' Run viper with decoupleR python script.
#'
#' To run decoupleR python version scanpy, pandas, and decoupleR python libraries need to be installed.
#'
#' @param object_path path to anndata object
#' @param regulon_path path to regulon table
#' @param out_path Output path to save results
#' @export
run_py_decoupler_viper <- function(object_path, regulon_path, out_path) {
  reticulate::source_python(system.file("python_scripts/run_decoupleR.py", package = "LR2TF"))
  run_decoupler_viper(object_path, regulon_path, out_path)
}

#' Run weighted mean with decoupleR python script.
#'
#' To run decoupleR python version scanpy, pandas, and decoupleR python libraries need to be installed.
#'
#' @param object_path path to anndata object
#' @param regulon_path path to regulon table
#' @param out_path Output path to save results
#' @export
run_py_decoupler_wmean <- function(object_path, regulon_path, out_path) {
  reticulate::source_python(system.file("python_scripts/run_decoupleR.py", package = "LR2TF"))
  run_decoupler_wmean(object_path, regulon_path, out_path)
}

#' Run ulm with decoupleR python script.
#'
#' To run decoupleR python version scanpy, pandas, and decoupleR python libraries need to be installed.
#'
#' @param object_path path to anndata object
#' @param regulon_path path to regulon table
#' @param out_path Output path to save results
#' @export
run_py_decoupler_ulm <- function(object_path, regulon_path, out_path) {
  reticulate::source_python(system.file("python_scripts/run_decoupleR.py", package = "LR2TF"))
  run_decoupler_ulm(object_path, regulon_path, out_path)
}

#' Run mlm with decoupleR python script.
#'
#' To run decoupleR python version scanpy, pandas, and decoupleR python libraries need to be installed.
#'
#' @param object_path path to anndata object
#' @param regulon_path path to regulon table
#' @param out_path Output path to save results
#' @export
run_py_decoupler_mlm <- function(object_path, regulon_path, out_path) {
  reticulate::source_python(system.file("python_scripts/run_decoupleR.py", package = "LR2TF"))
  run_decoupler_mlm(object_path, regulon_path, out_path)
}

#' Run mlm with decoupleR python script.
#'
#' To run decoupleR python version scanpy, pandas, and decoupleR python libraries need to be installed.
#'
#' @param object_path path to anndata object
#' @param condition_field name of meta-data field with condition information
#' @param cluster_field name of meta-data field with cluster information
#' @export
run_py_liana <- function(object_path, regulon_path, out_path) {
  reticulate::source_python(system.file("python_scripts/run_liana.py", package = "LR2TF"))
  run_liana(object_path, condition_field, cluster_field)
}

#' Loading old dorothea package regulons (default if no regulon is provided)
#'
#' @param arguments_list list of user defined arguments
#' @return val_arguments validated list of arguments
#' @export
validate_input_arguments <- function(arguments_list) {
  if (is.null(arguments_list$out_path)) {
    print("Please provide an output path")
  } else {
    if (substring(arguments_list$out_path, length(arguments_list$out_path) - 1, length(arguments_list$out_path)) != "/") {
      arguments_list$out_path <- paste0(arguments_list$out_path, "/")
    }
  }
  if (is.null(arguments_list$celltype)) {
    print("Please provide the name of the metadata field containing cell type annotations")
  }
  if (is.null(arguments_list$condition)) {
    print("Please provide the name of the metadata field containing condition annotations")
  }
  if (is.null(arguments_list$organism)) {
    arguments_list$organism <- "human"
  }
  if (is.null(arguments_list$comparison_list)) {
    arguments_list$comparison_list <- NA
  }
  if (is.null(arguments_list$logfc)) {
    arguments_list$logfc <- 0.0
  }
  if (is.null(arguments_list$pval)) {
    arguments_list$pval <- 0.05
  }
  if (is.null(arguments_list$reg)) {
    arguments_list$reg = load_dorothea_regulon(arguments_list$organism)
  } else {
     if (typeof(arguments_list$reg) == "character") {
       arguments_list$reg <-read.csv(arguments_list$reg, header = TRUE)
    }
    if(!all(c("source", "target", "weight") %in% names(arguments_list$reg))){
      stop("Not all necessary columns found in regulon table! Please make sure that the regulon has the columns source, target and weight!")
    }
  }
  return(arguments_list)
}

#' Run TF activity analysis
#'
#' Description
#'
#' @param seuratobject Input Seurat Object
#' @param tf_activities matrix with TF activities for each cell in the scRMA-seq data
#' @param arguments_list named R list with custom options for the analysis
#' @export
tf_activity_analysis <- function(seuratobject, tf_activities = NA, arguments_list) {

  if (typeof(seuratobject) == "character") {
    seuratobject <- readRDS(seuratobject)
  }

  arguments_list <- validate_input_arguments(arguments_list)

  Idents(object = seuratobject) <- arguments_list$celltype
  dir.create(arguments_list$out_path)
  tf_path <- paste0(arguments_list$out_path, "TF_results/")
  dir.create(tf_path)

  if (is.na(tf_activities)[[1]]) {
    seuratobject <- decoupleR_viper_analysis(seuratobject, arguments_list$reg)
  } else {
    if (typeof(tf_activities) == "character") {
      tf_activities <- t(read.csv(tf_activities, header = TRUE, row.names = 1))
    }

    tf_acticities <- CreateAssayObject(data = tf_activities)
    seuratobject[["tf_activities"]] <- tf_acticities
  }

  Idents(object = seuratobject) <- arguments_list$condition
  if (length(arguments_list$comparison_list) > 0 & length(levels(Idents(seuratobject))) < 2) {
    arguments_list$comparison_list <- NA
    print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")
  }
  Idents(object = seuratobject) <- arguments_list$celltype

  if (is.na(arguments_list$comparison_list)[[1]]) {
    seuratobject[['tf_annotation']] <- Idents(object = seuratobject)
    result_list <- list()
    gene_expression_list <- list()
    CTR_cluster_list <- list()
    intranet_cluster_list <- list()

    Idents(object = seuratobject) <- arguments_list$condition
    seuratobject_list <- SplitObject(seuratobject, split.by = arguments_list$condition)
    for (name in names(seuratobject_list)) {
      sub_object <- seuratobject_list[[name]]

      name <- str_replace_all(name, "[,;.:-]", "_")

      sub_object.averages <- AverageExpression(sub_object, group.by = arguments_list$celltype, assays = "RNA")
      # write.csv(seuratobject.averages[["RNA"]], file =
      #   paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
      #          name, '.csv'))

      tf_activity_scores = get_significant_tfs_single(sub_object, name, tf_path, pval = arguments_list$pval, log2fc = arguments_list$logfc)
      result_list[[name]] = tf_activity_scores
      gene_expression_list[[paste0(name, "_average_expression")]] = sub_object.averages[["RNA"]]

      if (arguments_list$organism == "human") {
        CTR_cluster_list[[name]] <- generate_CrossTalkeR_input(tf_activity_scores[["cluster"]],
                                                                                 gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                 arguments_list$reg)
      }else {
        CTR_cluster_list[[name]] <- generate_CrossTalkeR_input_mouse(tf_activity_scores[["cluster"]],
                                                                                       gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                       arguments_list$reg)
      }
      intranet_cluster_list[[name]] <- generate_intracellular_network(tf_activity_scores[["cluster"]],
                                                                      gene_expression_list[[paste0(name, "_average_expression")]],
                                                                      arguments_list$reg,
                                                                      arguments_list$organism)
    }
    tf <- new("TFObj",
              tf_activities_condition = list(),
              tf_activities_cluster = result_list,
              average_gene_expression = gene_expression_list,
              regulon = arguments_list$reg,
              CTR_input_condition = list(),
              CTR_input_cluster = CTR_cluster_list,
              intracellular_network_condition = list(),
              intracellular_network_cluster = intranet_cluster_list)

    saveRDS(tf, file = paste0(tf_path, "result_TF_object.RDS"))
    #saveRDS(seuratobject, file = paste0(tf_path, "TF_seurat_object.RDS"))
    return(tf)
  }

  else {
    out_path_compared <- paste0(tf_path, "compared")
    dir.create(out_path_compared)
    compared_significant_tfs <- condition_comparison_significant(seuratobject, out_path_compared, arguments_list$celltype, arguments_list$condition, arguments_list$comparison_list)

    plot_condition_tf_activities(compared_significant_tfs, out_path_compared)
    plot_condition_tf_activities_compressed(compared_significant_tfs, out_path_compared)

    seuratobject[['tf_annotation']] <- Idents(object = seuratobject)
    seuratobject_list <- SplitObject(seuratobject, split.by = arguments_list$condition)

    result_condition_list <- list()
    result_cluster_list <- list()
    gene_expression_list <- list()
    CTR_condition_list <- list()
    CTR_cluster_list <- list()
    intranet_condition_list <- list()
    intranet_cluster_list <- list()
    for (name in names(seuratobject_list)) {
      sub_object <- seuratobject_list[[name]]

      compared_tfs = data.frame(
        gene = character(),
        tag = character(),
        cluster = character()
      )

      for (result_name in names(compared_significant_tfs)) {
        if (grepl(name, result_name, fixed = TRUE)) {
          tf_condition_significant = compared_significant_tfs[[result_name]]
          #tf_condition_significant = tf_condition_significant[(tf_condition_significant$tag == "***"),]
          tf_condition_significant = filter(tf_condition_significant, FDR < as.double(arguments_list$pval))
          tf_condition_significant = filter(tf_condition_significant, logFC > as.double(arguments_list$logfc) | logFC < (0 - as.double(arguments_list$logfc)))
          tf_condition_significant = tf_condition_significant[c("tf", "tag", "CellType")]
          tf_condition_significant = tf_condition_significant %>% rename(gene = tf, cluster = CellType)
          compared_tfs = rbind(compared_tfs, tf_condition_significant)
        }
      }

      compared_tfs = compared_tfs[!duplicated(rownames(compared_tfs)),]

      name <- str_replace_all(name, "[,;.:-]", "_")

      sub_object.averages <- AverageExpression(sub_object, group.by = arguments_list$celltype, assays = "RNA")
      # write.csv(sub_object.averages[["RNA"]], file =
      #   paste0(arguments_list$out_path, 'average_gene_expression_by_cluster_',
      #          name, '.csv'))

      tf_activity_scores = get_significant_tfs(sub_object,
                                               name,
                                               tf_path,
                                               compared_tfs,
                                               pval = arguments_list$pval,
                                               log2fc = arguments_list$logfc)
      result_condition_list[[name]] <- tf_activity_scores[["condition"]]
      result_cluster_list[[name]] <- tf_activity_scores[["cluster"]]
      gene_expression_list[[paste0(name, "_average_expression")]] = sub_object.averages[["RNA"]]

      if (arguments_list$organism == "human") {
        CTR_condition_list[[name]] <- generate_CrossTalkeR_input(tf_activity_scores[["condition"]],
                                                                                   gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                   arguments_list$reg)
        CTR_cluster_list[[name]] <- generate_CrossTalkeR_input(tf_activity_scores[["cluster"]],
                                                                                 gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                 arguments_list$reg)
      }else {
        CTR_condition_list[[name]] <- generate_CrossTalkeR_input_mouse(tf_activity_scores[["condition"]],
                                                                                         gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                         arguments_list$reg)
        CTR_cluster_list[[name]] <- generate_CrossTalkeR_input_mouse(tf_activity_scores[["cluster"]],
                                                                                       gene_expression_list[[paste0(name, "_average_expression")]],
                                                                                       arguments_list$reg)
      }
      intranet_condition_list[[name]] <- generate_intracellular_network(tf_activity_scores[["condition"]],
                                                                        gene_expression_list[[paste0(name, "_average_expression")]],
                                                                        arguments_list$reg,
                                                                        arguments_list$organism)
      intranet_cluster_list[[name]] <- generate_intracellular_network(tf_activity_scores[["cluster"]],
                                                                      gene_expression_list[[paste0(name, "_average_expression")]],
                                                                      arguments_list$reg,
                                                                      arguments_list$organism)
    }

    tf <- new("TFObj",
              tf_activities_condition = result_condition_list,
              tf_activities_cluster = result_cluster_list,
              average_gene_expression = gene_expression_list,
              regulon = arguments_list$reg,
              CTR_input_condition = CTR_condition_list,
              CTR_input_cluster = CTR_cluster_list,
              intracellular_network_condition = intranet_condition_list,
              intracellular_network_cluster = intranet_cluster_list)

    saveRDS(tf, file = paste0(tf_path, "result_TF_object.RDS"))
    #saveRDS(seuratobject, file = paste0(tf_path, "TF_seurat_object.RDS"))

    return(tf)
  }

}
