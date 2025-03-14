
library(LR2TF)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(maditr)

TF_object <- readRDS("LR2TF_test_run/results/TF_results/result_TF_object.RDS")
TF_object@intracellular_network_condition
TF_object@intracellular_network_cluster
avg_gene_expr_ctrl <- TF_object@average_gene_expression$control_average_expression

head(avg_gene_expr_ctrl, 10)
avg_gene_expr_ctrl["IKZF1",]
write.csv(avg_gene_expr_ctrl, file = paste0('r_avg_expr_ctrl.csv'))

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
ligands_human <- load(file='ligands_human.rda')
str(ligands)


#################################


data <- readRDA("new_test//my_data.rda")


generate_CrossTalkeR_input <-
  function(tf_activities,
           gene_expression,
           regulon = NA,
           organism = "human") {
    if (any(tf_activities$z_score > 0)) {
      regulon <- regulon %>%
        rename(tf = source)

      organism <-"human"

      if(organism == "human") {
        ligands <- ligands_human
        R2TF <- aggregate(RTF_DB$receptor ~ RTF_DB$tf, FUN = c)
      } else if(organism == "mouse") {
        ligands <- ligands_mouse
        R2TF <- aggregate(RTF_DB_mouse$receptor ~ RTF_DB_mouse$tf, FUN = c)
      } else {
        stop("Invalid organism to generate CrossTalkeR input!")
      }

      colnames(R2TF) <- c('tf', 'receptors')
      R2TF <- R2TF %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')
      #print(R2TF)
      sorted_regulon <-
        aggregate(regulon$target ~ regulon$tf, FUN = c)
      colnames(sorted_regulon) <- c('tf', 'targets')
      sorted_regulon <- sorted_regulon %>%
        remove_rownames %>%
        tibble::column_to_rownames(var = 'tf')

      tf_activities <- tf_activities %>%
        filter(z_score > 0)

      output_df <- LR2TF::create_empty_CTR_dataframe()

      for (row in 1:nrow(tf_activities)) {

        r_tf <- LR2TF::create_empty_CTR_dataframe()
        tf_l <- LR2TF::create_empty_CTR_dataframe()

        #if (tf_activities[row, "z_score"] > 0) {
        tf <- as.character(tf_activities[row, "gene"])
        print(paste0("tf", tf))

        targets <- sorted_regulon[tf,][1]
        print(paste0("targets", targets[[1]]))
        #print(ligands)
        receptors <- R2TF[tf,][1]
        print(paste0("receptors", receptors))
        tf_ligands <- intersect(targets[[1]], ligands)
        print(paste0("tf_ligands", tf_ligands))
        if (length(tf_ligands) > 0) {
          for (ligand in tf_ligands) {
            expressed <- FALSE
            if (ligand %in% rownames(gene_expression)) {
              ex_value <- gene_expression[ligand, tf_activities[row, "cluster"]]
              if (ex_value != 0) {
                expressed <- TRUE
              }
            }

            if (expressed == TRUE) {
              df <- LR2TF::add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                               tf_activities[row, "cluster"],
                                               tf_activities[row, "gene"],
                                               ligand,
                                               'Transcription Factor',
                                               'Ligand',
                                               tf_activities[row, "z_score"])
              tf_l <- rbind(tf_l, df)
            }
          }
        }
        if (length(receptors[[1]]) > 0) {
          for (receptor in receptors[[1]]) {
            df <- LR2TF::add_entry_to_CTR_dataframe(tf_activities[row, "cluster"],
                                             tf_activities[row, "cluster"],
                                             receptor,
                                             tf_activities[row, "gene"],
                                             'Receptor',
                                             'Transcription Factor',
                                             tf_activities[row, "z_score"])
            r_tf <- rbind(r_tf, df)
          }
        }
        #}

        r_tf$gene_A <- gsub("_", "+", r_tf$gene_A, fixed = TRUE)
        r_tf$gene_B <- gsub("_", "+", r_tf$gene_B, fixed = TRUE)
        tf_l$gene_A <- gsub("_", "+", tf_l$gene_A, fixed = TRUE)
        tf_l$gene_B <- gsub("_", "+", tf_l$gene_B, fixed = TRUE)

        output_df <- rbind(output_df, r_tf)
        output_df <- rbind(output_df, tf_l)
      }
    } else {
      output_df = NA
    }

    return(output_df)
    print(output_df)
  }


library(dplyr)
library(tibble)

#tf_activities <- results@tf_activities_condition$control
tf_activities <- read.csv("script_test/TF_results/control/significant_condition_tf_results_control.csv")
tf_activities <- tf_activities %>% rename_at('t_value', ~'z_score')
tf_activities <- tf_activities[tf_activities$cluster == "Fibroblast",]
print(tf_activities)


#tf_activities <- read.csv("script_test/TF_results/PMF_MF2/significant_condition_tf_results_PMF,MF2.csv")
#tf_activities <- tf_activities %>% rename_at('t_value', ~'z_score')
print(tf_activities)

gene_expression <- results@average_gene_expression$control_average_expression
print(gene_expression)

regulon <- read.csv("LR2TF_test_run/filterd_regulon.csv")
names(regulon)[names(regulon) == 'source'] <- 'tf'
colnames(regulon) <- c('tf', 'target')

RTF_DB <- read.csv("rtf_db_human.csv")


ligands_human <- read.csv("ligands_human.csv")


result_r <- generate_CrossTalkeR_input(tf_activities, gene_expression, regulon, organism = "human")
print(result_r)

write.csv(result_r, "r_output.csv", row.names = FALSE)
write.table(result_r, file = "r_output.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

result_py <- generate_CrossTalkeR_input(tf_activities, gene_expression, regulon, organism = "human")
write.table(result_py, file = "py_output_control_cond_in_R.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
print(result_py)

table_ctr <- read.csv("LR2TF_test_run/CTR_LR.csv")
table_exp <- read.csv("/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/EXP_LR.csv")

ctr_input_py <- LR2TF::combine_LR_and_TF_complexes(result_py, table_ctr, "z", "control_cond_py_in_R_try2")
exp_input_py <- LR2TF::combine_LR_and_TF_complexes(result_py, table_exp, "z", "PMF_cond_py_in_R")

load(file="RTF_DB_2.rda")
head(rtf_db_R[1])
print(RTF_DB_2)
RTF_DB <- RTF_DB_2
write.csv(RTF_DB_2, "rtf_db_human.csv", row.names = TRUE)

load("ligands_human.rda")
write.csv(ligands_human, "ligands_human_diff_maybe.csv", row.names = TRUE)


result_py_py <- read.csv("py_ctr_input_wo_ctr_exp_tables.csv")
ctr_input_py <- LR2TF::combine_LR_and_TF_complexes(result_py_py, table_ctr, "z", "contr_cond_CTR_input_py_but_combine_LR_TF_with_R")
