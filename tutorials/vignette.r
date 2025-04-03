library(LR2TF)
#library(SeuratObject)
library(Seurat)

#remotes::install_github("CostaLab/LR2TF", ref="dev_cleanup", force =TRUE)


#test dataset from package:
data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example


parameters <- list("out_path" = "/home/larissa/Documents/LR2TF_HiWi/new_test/",
                   reg = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/filterd_regulon.csv",
                   "organism" = "human",
                   "celltype" = "new_annotation", #name of the meta data field defining cell identities
                   "condition" = "protocol", #name of the meta data field defining conditions
                   "comparison_list" = list(c("PMF,MF2", "control")), #list of condition comparison to consider
                   "logfc" = 0.5,
                   "pval" = 0.05) #thresholds for logfc and pval used in differential transcription factor analysis



results <- LR2TF::IntraTalker_analysis(seuratobject = seurat_object,
                                       tf_activities = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/decoupler_results.csv",
                                       arguments_list = parameters)

#1. tf_activities_condition -> tables with condition significant transcription factors for each compared condition
#2. tf_activities_cluster -> tables with cluster specific transcription factors for all conditions in the data
#3. average_gene_expression -> matrices for each condition with average gene expressions
#4. regulon -> regulon used for the analysis as specified by the user
#5. CTR_input_condition -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on condition specific transcription factors (input for CrossTalker)
#6. CTR_input_cluster -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on cluster specific transcription factors (input for CrossTalker)
#7. intracellular_network_condition -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on condition specific transcription factors
#8. intracellular_network_cluster -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on cluster specific transcription factors


write.csv(results@CTR_input_condition[["control"]], "R_ctr_input_wo_exp_ctr_tables.csv")

write.csv(results@intracellular_network_condition[["control"]], "R_intra_network_ctrl.csv")
write.csv(results@intracellular_network_condition[["PMF_MF2"]], "R_intra_network_PMF.csv")

write.csv(results@intracellular_network_cluster[["control"]], "R_intra_network_ctrl_cluster.csv")
write.csv(results@intracellular_network_cluster[["PMF_MF2"]], "R_intra_network_PMF_cluster.csv")

table_ctr <- read.csv("/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/CTR_LR.csv", row.names = NULL)
table_exp <- read.csv("/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/EXP_LR.csv")

table_ctr$X <- NULL
table_exp$X <- NULL

ctr_inptu <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")


ctr_input_cluster <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_cluster[["control"]], table_ctr, parameters$out_path, "control_cluster")
exp_input_cluster <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_cluster[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2_cluster")

########################
#WINDOWS:

library(LR2TF)
library(Seurat)

#test dataset from package:
data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example


parameters <- list("out_path" = "new_test\\",
                   "reg" = "LR2TF_test_run\\filterd_regulon.csv",
                   "organism" = "human",
                   "celltype" = "new_annotation", #name of the meta data field defining cell identities
                   "condition" = "protocol", #name of the meta data field defining conditions
                   "comparison_list" = list(c("PMF,MF2", "control")), #list of condition comparison to consider
                   "logfc" = 0.5,
                   "pval" = 0.05) #thresholds for logfc and pval used in differential transcription factor analysis



results <- LR2TF::IntraTalker_analysis(seuratobject = seurat_object,
                                       tf_activities = "LR2TF_test_run\\decoupler_results.csv",
                                       arguments_list = parameters)


#1. tf_activities_condition -> tables with condition significant transcription factors for each compared condition
#2. tf_activities_cluster -> tables with cluster specific transcription factors for all conditions in the data
#3. average_gene_expression -> matrices for each condition with average gene expressions
#4. regulon -> regulon used for the analysis as specified by the user
#5. CTR_input_condition -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on condition specific transcription factors (input for CrossTalker)
#6. CTR_input_cluster -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on cluster specific transcription factors (input for CrossTalker)
#7. intracellular_network_condition -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on condition specific transcription factors
#8. intracellular_network_cluster -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on cluster specific transcription factors

save(results, file='new_test\\my_data.rda')
results <- readRDA(file='new_test\\my_data.rda')

write.csv(results@CTR_input_condition[["control"]], "R_ctr_input_wo_exp_ctr_tables.csv")

#print(results@intracellular_network_condition[["control"]])
write.csv(results@intracellular_network_condition[["control"]], "R_intra_network_ctrl.csv")
write.csv(results@intracellular_network_condition[["PMF_MF2"]], "R_intra_network_PMF.csv")

write.csv(results@intracellular_network_cluster[["control"]], "R_intra_network_ctrl_cluster.csv")
write.csv(results@intracellular_network_cluster[["PMF_MF2"]], "R_intra_network_PMF_cluster.csv")

table_ctr <- read.csv("LR2TF_test_run\\CTR_LR.csv")
table_exp <- read.csv("LR2TF_test_run\\EXP_LR.csv")


table_ctr$X <- NULL
table_exp$X <- NULL
ctr_inptu <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")

ctr_input_cluster <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_cluster[["control"]], table_ctr, parameters$out_path, "control_cluster")
exp_input_cluster <- LR2TF::combine_LR_and_TF_complexes(results@CTR_input_cluster[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2_cluster")



####################################################
library(dplyr)
library(tibble)
library(tidyr)

result_py <- read.csv("py_output_in_R.csv")

#condition
  #tf_table <- results@CTR_input_condition[["control"]]
  tf_table <- result_py
  lr_table <- table_ctr
  intra_connections <- tf_table[NULL, ]
  for (celltype in unique(append(lr_table$source, lr_table$target))) {
    lr_filtered_ligands <- lr_table[lr_table$source == celltype, ]
    lr_filtered_receptors <- lr_table[lr_table$target == celltype, ]
    lr_ligands <- unique(lr_filtered_ligands$gene_A)
    lr_receptors <- unique(lr_filtered_receptors$gene_B)
    print(celltype)
    print(lr_receptors)
    tf_table_receptors <- tf_table[tf_table$target == celltype & tf_table$type_gene_A == "Receptor", ]
    tf_table_ligands <- tf_table[tf_table$source == celltype & tf_table$type_gene_B == "Ligand", ]
    tf_receptor_interactions <- tf_table_receptors %>%
      filter(gene_A %in% lr_receptors)
    tf_ligand_interactions <- tf_table_ligands %>%
      filter(gene_B %in% lr_ligands)
    intra_connections <- rbind(intra_connections, tf_receptor_interactions, tf_ligand_interactions)
  }
  intra_connections$all_pair <- paste0(
    intra_connections$source, "/",
    intra_connections$gene_A, "/",
    intra_connections$target, "/",
    intra_connections$gene_B
  )
  intra_connections <- intra_connections[!duplicated(intra_connections$all_pair), ]
  
  intra_connections$all_pair <- NULL
  complete_interactions <- rbind(intra_connections, lr_table)
  
  write.csv(complete_interactions, paste0("new_test/CrossTalkeR_input_control.csv"), row.names = FALSE)
  write.csv(complete_interactions, paste0("CrossTalkeR_input_control_py_in_R.csv"), row.names = FALSE)

#cluster


  tf_table <- results@CTR_input_cluster[["control"]]
  #print(tf_table)
  lr_table <- table_ctr
  intra_connections <- tf_table[NULL, ]
  for (celltype in unique(append(lr_table$source, lr_table$target))) {
    lr_filtered_ligands <- lr_table[lr_table$source == celltype, ]
    lr_filtered_receptors <- lr_table[lr_table$target == celltype, ]
    lr_ligands <- unique(lr_filtered_ligands$gene_A)
    lr_receptors <- unique(lr_filtered_receptors$gene_B)
    print(celltype)
    print(lr_receptors)
    tf_table_receptors <- tf_table[tf_table$target == celltype & tf_table$type_gene_A == "Receptor", ]
    tf_table_ligands <- tf_table[tf_table$source == celltype & tf_table$type_gene_B == "Ligand", ]
    tf_receptor_interactions <- tf_table_receptors %>%
      filter(gene_A %in% lr_receptors)
    tf_ligand_interactions <- tf_table_ligands %>%
      filter(gene_B %in% lr_ligands)
    intra_connections <- rbind(intra_connections, tf_receptor_interactions, tf_ligand_interactions)
  }
  intra_connections$all_pair <- paste0(
    intra_connections$source, "/",
    intra_connections$gene_A, "/",
    intra_connections$target, "/",
    intra_connections$gene_B
  )
  intra_connections <- intra_connections[!duplicated(intra_connections$all_pair), ]
  
  intra_connections$all_pair <- NULL
  complete_interactions <- rbind(intra_connections, lr_table)
  
  write.csv(complete_interactions, paste0("new_test/CrossTalkeR_input_control_cluster.csv"), row.names = FALSE)
  write.csv(complete_interactions, paste0("CrossTalkeR_input_control_cluster_py_in_R.csv"), row.names = FALSE)
