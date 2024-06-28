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
data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example

head(seurat_object@assays$RNA@data)
seurat_object
#sort(seurat_object@assays$RNA@data, decreasing = TRUE)


#tmp=seurat_object@assays$RNA@counts
#seurat_object@assays$RNA@counts=seurat_object@assays$RNA@data

#sub_object.averages <- AverageExpression(seurat_object, group.by = "new_annotation", assay = "RNA", slot = "count") # nolint: line_length_linter.
#as.data.frame(sub_object.averages)

#write.csv(sub_object.averages["RNA"], file = "R_avg_2.csv")
#write.csv(seurat_object.averages["RNA"], file = "R_avg.csv")

#LR2TF::convert_seurat_to_anndata(seurat_object, "D:\\studium\\HiWi\\LR2TF\\")


#The file "anndata_object.h5ad" will be saved into the user defined path and can then be used to perform the predictions. Beside the scRNA-seq data file, we also need to define a regulon database in form of a csv file with the coloumns "source", "target" and "weight". Within this package we provide the dorothea databases for human and mouse, downloaded from the decoupleR package. These files also contain the column "confidence" (levels A to D) with information on how well described a transcription factor and target gene interaction is in different resources. We recommend using the confidence levels A and B.



regulon <- read.csv(system.file("regulons", "human_dorothea_reg.csv", package = "LR2TF"), row.names = 1)
filtered_regulon <- regulon[regulon$confidence %in% c("A","B"),]
write.csv(regulon, "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/filterd_regulon.csv")

######

parameters <- list("out_path" = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/results/",
                   reg = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/filterd_regulon.csv",
                   "organism" = "human",
                   "celltype" = "new_annotation", #name of the meta data field defining cell identities
                   "condition" = "protocol", #name of the meta data field defining conditions
                   "comparison_list" = list(c("PMF,MF2", "control")) , #list of condition comparison to consider
                   "logfc" = 0.5,
                   "pval" = 0.05) #thresholds for logfc and pval used in differential transcription factor analysis


results <- LR2TF::tf_activity_analysis(seuratobject = seurat_object,
                                       tf_activities = "NA",
                                       arguments_list = parameters)




#table_ctr <- read.csv("D:\\studium\\HiWi\\LR2TF\\CTR_LR.csv", row.names = 1)
#table_exp <- read.csv("D:\\studium\\HiWi\\LR2TF\\EXP_LR.csv", row.names = 1)

#ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
#exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")

ctr_file <- "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/CTR_LR.csv"
exp_file <- "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/EXP_LR.csv"

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], ctr_file, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], exp_file, parameters$out_path, "PMF_MF2")

write.csv(ctr_inptu,"/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/ctr_lr_tf.csv", row.names = FALSE)
write.csv(exp_input,"/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_single_cond/exp_lr_tf.csv", row.names = FALSE)