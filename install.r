if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
BiocManager::install("OmnipathR")
BiocManager::install("dorothea")
BiocManager::install("viper")

remotes::install_github("cellgeni/sceasy")
remotes::install_github("CostaLab/LR2TF", build_vignettes = FALSE)
remotes::install_github("satijalab/seurat")
install.packages("anndata")

library(LRTF2)
library(Seurat)
library(anndata)

data(bone_marrow_stromal_cell_example, package = "LR2TF")
seurat_object <- bone_marrow_stromal_cell_example

#View(seurat_object)


LR2TF::convert_seurat_to_anndata(seurat_object, "D:\\studium\\HiWi\\")
```

#The file "anndata_object.h5ad" will be saved into the user defined path and can then be used to perform the predictions. Beside the scRNA-seq data file, we also need to define a regulon database in form of a csv file with the coloumns "source", "target" and "weight". Within this package we provide the dorothea databases for human and mouse, downloaded from the decoupleR package. These files also contain the column "confidence" (levels A to D) with information on how well described a transcription factor and target gene interaction is in different resources. We recommend using the confidence levels A and B.

```{r, include = TRUE, eval = FALSE}

regulon <- read.csv(system.file("regulons", "human_dorothea_reg.csv", package = "LR2TF"), row.names = 1)
filtered_regulon <- regulon[regulon$confidence %in% c("A","B"),]
write.csv(regulon, "filterd_regulon.csv")
