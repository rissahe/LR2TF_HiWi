
#install.packages("jsonlite")
#library(jsonlite)
#library(BiocManager)
#install.packages("fgsea")
#BiocManager::install(version = "3.20")
#BiocManager::install("fgsea")
#BiocManager::install("clusterProfiler")
#library(clusterProfiler)
pandoc_install("default")
pandoc_version()

library(devtools)
install_github("CostaLab/CrossTalkeR", ref="dev")
remove.packages("BiocManager")

library(CrossTalkeR)
library(igraph)
library(stringr)
library(tibble)
#library(pandoc)

ctr_input <- read.csv("/home/larissa/Documents/LR2TF_HiWi/script_test/CrossTalkeR_input_control.csv")
exp_input <- read.csv("/home/larissa/Documents/LR2TF_HiWi/script_test/CrossTalkeR_input_PMF_MF2.csv")

ctr_input$Unnamed..0<- NULL
exp_input$Unnamed..0<- NULL


paths <- c(
  "control" = "script_test/CrossTalkeR_input_control.csv",
  "PMF_MF2" = "script_test/CrossTalkeR_input_PMF_MF2.csv"
)


conds <- names(paths)
is.character(paths[[conds[1]]])
#yes was considered character even when reading the csvs beforehand amd r showing the objects as dataframe

output <- ("ctr_results")
data <- generate_report(paths,
  out_path = paste0(output, "/"),
  out_file = "HumanMyfib_TF_example.html",
  output_fmt = "html_document",
  report = TRUE,
  comparison = list(c("PMF_MF2", "control"))
)


plot_bar_rankings(data, "PMF_MF2_x_control_ggi", ranking = "Mediator", type = "TF")
dev.off()

pagerank_table = data@rankings$PMF_MF2_x_control_ggi %>%
  select(c("nodes", "Pagerank"))
pagerank_table = as.data.frame(pagerank_table)
rownames(pagerank_table) = pagerank_table$nodes

plot_graph_sankey_tf(data@tables$PMF_MF2_x_control,
    pagerank_table,
    target = "TGFB1|L",
    cluster = "Megakaryocyte",
    target_type = "Ligand",
    plt_name = "Signaling Upstream TGFB1 in Megakaryocytes")
dev.off()

LR_workspace <- readRDS("ctr_results/LR_data_final.Rds")
LR_workspace@tables$pagerank_table
