library(LR2TF)
library(Seurat)


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



results <- LR2TF::tf_activity_analysis(seuratobject = seurat_object,
                                       tf_activities = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/decoupler_results.csv",
                                       arguments_list = parameters)

1. tf_activities_condition -> tables with condition significant transcription factors for each compared condition
2. tf_activities_cluster -> tables with cluster specific transcription factors for all conditions in the data
3. average_gene_expression -> matrices for each condition with average gene expressions
4. regulon -> regulon used for the analysis as specified by the user
5. CTR_input_condition -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on condition specific transcription factors (input for CrossTalker)
6. CTR_input_cluster -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on cluster specific transcription factors (input for CrossTalker)
7. intracellular_network_condition -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on condition specific transcription factors
8. intracellular_network_cluster -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on cluster specific transcription factors


table_ctr <- read.csv("/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/CTR_LR.csv")
table_exp <- read.csv("/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/EXP_LR.csv")

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")

ctr_file <- "/home/larissa/Documents/Larissa_HiWi/LR2TF_test/control_lr_results.csv"
exp_file <- "/home/larissa/Documents/Larissa_HiWi/LR2TF_test/PMF,MF2_lr_results.csv"

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], ctr_file, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], exp_file, parameters$out_path, "PMF_MF2")

```
