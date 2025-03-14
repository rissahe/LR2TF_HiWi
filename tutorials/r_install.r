#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#BiocManager::install("OmnipathR")
#BiocManager::install("dorothea")
#BiocManager::install("viper")
#BiocManager::install("scran")

#install.packages("devtools")
#install.packages("rcompanion")
install.packages("spatstat.utils")
#install.packages("anndata")

remove.packages("spatstat.utils")
#anndata::install_anndata()
#remotes::install_github("cellgeni/sceasy")

#pandoc/VignetteBuilder is necessary for vignettes. install knitr?
# install.packages("pandoc") is not correct


install.packages("xfun")

remotes::install_github("CostaLab/LR2TF", build_vignettes = FALSE)
   
library(Seurat)
library(LR2TF)
library(Exact)
---
title: "LR2TF-Usage"
output: html_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, include = TRUE, eval = FALSE}
library(LR2TF)
library(Seurat)

#load own dataset:
#seurat_object <- readRDS("/path/to/seurat_object.Rds")

#test dataset from package:
data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example

LR2TF::convert_seurat_to_anndata(seurat_object, "/home/larissa/Documents/Larissa_HiWi/LR2TF_test/")
```

The file "anndata_object.h5ad" will be saved into the user defined path and can then be used to perform the predictions. Beside the scRNA-seq data file, we also need to define a regulon database in form of a csv file with the coloumns "source", "target" and "weight". Within this package we provide the dorothea databases for human and mouse, downloaded from the decoupleR package. These files also contain the column "confidence" (levels A to D) with information on how well described a transcription factor and target gene interaction is in different resources. We recommend using the confidence levels A and B.

```{r, include = TRUE, eval = FALSE}

regulon <- read.csv(system.file("regulons", "human_dorothea_reg.csv", package = "LR2TF"), row.names = 1)
filtered_regulon <- regulon[regulon$confidence %in% c("A","B"),]
write.csv(filtered_regulon, "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/filterd_regulon.csv")

```

Now it is possible to run the transcription factor activity prediction with decoupleR and the uml method:
```{python, include = TRUE, eval = FALSE}


import scanpy as sc
import decoupler as dc
import pandas as pd

ann_data = sc.read_h5ad("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/anndata_object.h5ad")
reg = pd.read_csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/filterd_regulon.csv")

dc.run_ulm( mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

estimates =ann_data.obs['ulm_estimate']
estimates.to_csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/decoupler_results.csv")

```

## Using the LR2TF package
Now, that we have a transcription factor activity matrix, we can continue the analysis with the LR2TF package. In this case we will use our test data for this example. First of all, it is necessary to define the following parameters in form of a list:

```{r eval = FALSE}
parameters <- list("out_path" = "/home/larissa/Documents/LR2TF_HiWi/new_test/",
                   "reg" = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/filterd_regulon.csv",
                   "organism" = "human",
                   "celltype" = "new_annotation", 
                   "condition" = "protocol", 
                   "comparison_list" = list(c("PMF,MF2", "control")), 
                   "logfc" = 0.0,
                   "pval" = 0.05) 
```

After defining the necessary parameter the transcription factor activity can be performed by calling:

```{r, include = TRUE, eval = FALSE}
results <- LR2TF::tf_activity_analysis(seuratobject = seurat_object,
                                       tf_activities = "/home/larissa/Documents/LR2TF_HiWi/LR2TF_test_run/decoupler_results.csv",
                                       arguments_list = parameters)
```

The "results" object contains the results of the performed analyses, consisting of multiple tables inside the object:

1. tf_activities_condition -> tables with condition significant transcription factors for each compared condition
2. tf_activities_cluster -> tables with cluster specific transcription factors for all conditions in the data
3. average_gene_expression -> matrices for each condition with average gene expressions
4. regulon -> regulon used for the analysis as specified by the user
5. CTR_input_condition -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on condition specific transcription factors (input for CrossTalker)
6. CTR_input_cluster -> for each condition a table with receptor-transcription factor and transcription factor-ligand interactions based on cluster specific transcription factors (input for CrossTalker)
7. intracellular_network_condition -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on condition specific transcription factors
8. intracellular_network_cluster -> for each condition a table with receptor-transcription factor and transcription factor-target gene interactions based on cluster specific transcription factors

(Note that special characters might be exchanged by underscores, if they cause problems with the naming of the tables; eg PMF,MF2 -> PMF_MF2)

The last step is to combine previous results from ligand receptor interaction analyses (e.g. CellPhoneDB) with the transcription factor results. (In the case of the test data the ligand-receptor interactions are provided within the CrossTalkeR package.)
```{r, include = TRUE, eval = FALSE}
table_ctr <- read.csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test_run/CTR_LR.csv", row.names = 1)
table_exp <- read.csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test_run/EXP_LR.csv", row.names = 1)

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], table_ctr, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], table_exp, parameters$out_path, "PMF_MF2")

ctr_file <- "/home/larissa/Documents/Larissa_HiWi/LR2TF_test/control_lr_results.csv"
exp_file <- "/home/larissa/Documents/Larissa_HiWi/LR2TF_test/PMF,MF2_lr_results.csv"

ctr_inptu <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["control"]], ctr_file, parameters$out_path, "control")
exp_input <- LR2TF::combine_LR_and_TF(results@CTR_input_condition[["PMF_MF2"]], exp_file, parameters$out_path, "PMF_MF2")
```



