if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
BiocManager::install("OmnipathR")
BiocManager::install("dorothea")
BiocManager::install("viper")
BiocManager::install("scran")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

remotes::install_github("cellgeni/sceasy")
remotes::install_github("CostaLab/LR2TF", build_vignettes = FALSE)
remotes::install_github("satijalab/seurat")
remotes::install_github("CostaLab/CrossTalkeR", build_vignettes = FALSE)
install.packages("anndata")
install.packages("rcompanion")
install.packages("DescTools")
install.packages("TH.data")

install.packages("devtools")
devtools::install_github('immunogenomics/presto')

library(LR2TF)
library(Seurat)
library(anndata)
library(devtools)
data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example


LR2TF::convert_seurat_to_anndata(seurat_object, "D:\\studium\\HiWi\\LR2TF\\")


#The file "anndata_object.h5ad" will be saved into the user defined path and can then be used to perform the predictions. Beside the scRNA-seq data file, we also need to define a regulon database in form of a csv file with the coloumns "source", "target" and "weight". Within this package we provide the dorothea databases for human and mouse, downloaded from the decoupleR package. These files also contain the column "confidence" (levels A to D) with information on how well described a transcription factor and target gene interaction is in different resources. We recommend using the confidence levels A and B.



regulon <- read.csv(system.file("regulons", "human_dorothea_reg.csv", package = "LR2TF"), row.names = 1)
filtered_regulon <- regulon[regulon$confidence %in% c("A","B"),]
write.csv(regulon, "LR2TF\\filterd_regulon.csv")

######

parameters <- list("out_path" = "D:\\studium\\HiWi\\LR2TF\\results",
                   reg = "D:\\studium\\HiWi\\LR2TF\\filterd_regulon.csv",
                   "organism" = "human",
                   "celltype" = "new_annotation", #name of the meta data field defining cell identities
                   "condition" = "protocol", #name of the meta data field defining conditions
                   "comparison_list" = list(c("PMF,MF2", "control")), #list of condition comparison to consider
                   "logfc" = 0.5,
                   "pval" = 0.05) #thresholds for logfc and pval used in differential transcription factor analysis


results <- LR2TF::tf_activity_analysis(seuratobject = seurat_object,
                                       tf_activities = "D:\\studium\\HiWi\\LR2TF\\decoupler_results.csv",
                                       arguments_list = parameters)




#table_ctr <- read.csv("D:\\studium\\HiWi\\LR2TF\\CTR_LR.csv", row.names = 1)
#table_exp <- read.csv("D:\\studium\\HiWi\\LR2TF\\EXP_LR.csv", row.names = 1)

#ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
#exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")

ctr_file <- "D:\\studium\\HiWi\\LR2TF\\CTR_LR.csv"
exp_file <- "D:\\studium\\HiWi\\LR2TF\\EXP_LR.csv"

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], ctr_file, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], exp_file, parameters$out_path, "PMF_MF2")

write.csv(ctr_inptu,"D:\\studium\\HiWi\\LR2TF\\ctr_lr_tf.csv", row.names = FALSE)
write.csv(exp_input,"D:\\studium\\HiWi\\LR2TF\\exp_lr_tf.csv", row.names = FALSE)