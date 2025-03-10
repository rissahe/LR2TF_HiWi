#install.packages('VennDiagram')
library(VennDiagram)

#########
#CTRL
#########

csv1 <- read.csv("new_test/CrossTalkeR_input_control.csv")
#row.names(csv1) <- NULL
csv1 <- csv1[c(1:193),]
csv1 <- csv1[csv1$MeanLR > 0,]


#csv1 <- read.csv("R_ctr_input_wo_exp_ctr_tables.csv")
#csv1 <- results@CTR_input_condition[["control"]]

csv1 <- csv1[c("source", "gene_A", "gene_B")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}

head(csv1_list)
SET1 <- csv1_list
length(csv1_list)

csv2 <- read.csv("script_test/CrossTalkeR_input_control.csv")
row.names(csv2) <- NULL
csv2 <- csv2[c(1:342),]
csv2 <- csv2[csv2$MeanLR > 0,]

#csv2 <- read.csv("py_ctr_input_wo_ctr_exp_tables.csv")

csv2 <- csv2[c("source", "gene_A", "gene_B")]


csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}

head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#genes_filtered_w_LR_table_no_neg_score & R results main scaled too
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_unique_CTR_control_condition.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/Py_unique_CTR_control_condition.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

#filtered_with_LR_table_no_neg_score_
pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_CTR_control_condition.pdf")
grid.draw(v)
dev.off()

################

#PMF
####################

csv1 <- read.csv("new_test/CrossTalkeR_input_PMF_MF2.csv")
row.names(csv1) <- NULL
csv1 <- csv1[c(1:269),]
csv1 <- csv1[csv1$MeanLR > 0,]

#csv1 <- read.csv("R_ctr_input_wo_exp_ctr_tables.csv")
#csv1 <- results@CTR_input_condition[["control"]]

csv1 <- csv1[c("source", "gene_A", "gene_B")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}

head(csv1_list)
SET1 <- csv1_list


csv2 <- read.csv("script_test/CrossTalkeR_input_PMF_MF2.csv")
row.names(csv2) <- NULL
csv2 <- csv2[c(1:277),]
csv2 <- csv2[csv2$MeanLR > 0,]

#csv2 <- read.csv("py_ctr_input_wo_ctr_exp_tables.csv")

csv2 <- csv2[c("source", "gene_A", "gene_B")]


csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}

head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_unique_CTR_PMF_condition.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/Py_unique_CTR_PMF_condition.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_CTR_input_PMF_condition.pdf")
grid.draw(v)
dev.off()

########################################################
#CRT input cluster control
#R and py are both scaled in main function

csv1 <- read.csv("new_test/CrossTalkeR_input_control_cluster.csv")
#row.names(csv1) <- NULL
csv1 <- csv1[c(1:754),]
csv1 <- csv1[csv1$MeanLR > 0,]


#csv1 <- read.csv("R_ctr_input_wo_exp_ctr_tables.csv")
#csv1 <- results@CTR_input_condition[["control"]]

csv1 <- csv1[c("source", "gene_A", "gene_B")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}

head(csv1_list)
SET1 <- csv1_list
length(csv1_list)

csv2 <- read.csv("script_test/CrossTalkeR_input_control_cluster.csv")
row.names(csv2) <- NULL
csv2 <- csv2[c(1:1097),]
csv2 <- csv2[csv2$MeanLR > 0,]

#csv2 <- read.csv("py_ctr_input_wo_ctr_exp_tables.csv")

csv2 <- csv2[c("source", "gene_A", "gene_B")]


csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}

head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#genes_filtered_w_LR_table_no_neg_score
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_unique_CTR_control_cluster.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/Py_unique_CTR_control_cluster.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

#filtered_with_LR_table_no_neg_score_
pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_CTR_input_control_cluster.pdf")
grid.draw(v)
dev.off()


########################################################
#CRT input cluster PMF MF2
#R and py are both scaled in main function

csv1 <- read.csv("new_test/CrossTalkeR_input_PMF_MF2_cluster.csv")
#row.names(csv1) <- NULL
csv1 <- csv1[c(1:741),]
csv1 <- csv1[csv1$MeanLR > 0,]


#csv1 <- read.csv("R_ctr_input_wo_exp_ctr_tables.csv")
#csv1 <- results@CTR_input_condition[["control"]]

csv1 <- csv1[c("source", "gene_A", "gene_B")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}

head(csv1_list)
SET1 <- csv1_list
length(csv1_list)

csv2 <- read.csv("script_test/CrossTalkeR_input_PMF_MF2_cluster.csv")
row.names(csv2) <- NULL
csv2 <- csv2[c(1:1747),]
csv2 <- csv2[csv2$MeanLR > 0,]

#csv2 <- read.csv("py_ctr_input_wo_ctr_exp_tables.csv")

csv2 <- csv2[c("source", "gene_A", "gene_B")]


csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}

head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#genes_filtered_w_LR_table_no_neg_score
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_unique_PMF_cluster.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/Py_unique_PMF_cluster.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

#filtered_with_LR_table_no_neg_score_
pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_CTR_input_PMF_cluster.pdf")
grid.draw(v)
dev.off()


##############################################################
#type gene = TF
#CTRL
################

csv1 <- read.csv("new_test/TF_results/control/significant_condition_tf_results_control.csv")
colnames(csv1)
csv1 <- csv1[csv1$z_score > 0,]

csv1 <- csv1[c("gene", "cluster")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}
length(csv1_list)
head(csv1_list)
SET1 <- csv1_list

csv2 <- read.csv("script_test/TF_results/control/significant_condition_tf_results_control.csv")

csv2 <- csv2[c("gene", "cluster")]

csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}
length(csv2_list)
head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#also no negative z score and filtered via LR table
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_significant_TF_control_condition.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_significant_TF_control_condition.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_significant_TF_control_condition.pdf")
grid.draw(v)
dev.off()
########################################################
#type =TF
#PMF.MF2
#####################

csv1 <- read.csv("new_test/TF_results/PMF_MF2/significant_condition_tf_results_PMF_MF2.csv")
csv1 <- csv1[csv1$z_score > 0,]
csv1 <- csv1[c("gene", "cluster")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}
length(csv1_list)
head(csv1_list)
SET1 <- csv1_list


csv2 <- read.csv("script_test/TF_results/PMF_MF2/significant_condition_tf_results_PMF,MF2.csv")

csv2 <- csv2[c("gene", "cluster")]

csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}
length(csv2_list)
head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#also no negative z score and filtered via LR table
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_significant_TF_PMF_condition.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_significant_TF_PMF_condition.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_significant_TF_PMF_condition.pdf")
grid.draw(v)
dev.off()

#####################################################################################
#tfs cluster 
#CTRL

csv1 <- read.csv("new_test/TF_results/control/significant_cluster_tf_results_control.csv")
colnames(csv1)
csv1 <- csv1[csv1$z_score > 0,]

csv1 <- csv1[c("gene", "cluster")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}
length(csv1_list)
head(csv1_list)
SET1 <- csv1_list

csv2 <- read.csv("script_test/TF_results/control/significant_cluster_tf_results_wilcoxon_control.csv")

csv2 <- csv2[c("gene", "cluster")]

csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}
length(csv2_list)
head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#also no negative z score and filtered via LR table
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_significant_TF_control_cluster.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_significant_TF_control_cluster.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_significant_TF_control_cluster.pdf")
grid.draw(v)
dev.off()
########################################################
#type =TF cluster
#PMF.MF2
#####################

csv1 <- read.csv("new_test/TF_results/control/significant_cluster_tf_results_control.csv")
csv1 <- csv1[csv1$z_score > 0,]
csv1 <- csv1[c("gene", "cluster")]

csv1_list <- list()
for (i in 1:nrow(csv1)) {
    row <- paste(csv1[i, ], collapse = ",")
    csv1_list <- append(csv1_list, row)
}
length(csv1_list)
head(csv1_list)
SET1 <- csv1_list


csv2 <- read.csv("script_test/TF_results/PMF_MF2/significant_cluster_tf_results_wilcoxon_PMF_MF2.csv")

csv2 <- csv2[c("gene", "cluster")]

csv2_list <- list()
for (i in 1:nrow(csv2)) {
    row <- paste(csv2[i, ], collapse = ",")
    csv2_list <- append(csv2_list, row)
}
length(csv2_list)
head(csv2_list)
SET2 <- csv2_list

v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R" , "Set PY "),
  filename = NULL,
  output = TRUE)


grid.newpage()
grid.draw(v)


#v[[5]]$label  <- paste(setdiff(SET1, intersect(SET1,SET2)), collapse="\n") 
#v[[6]]$label <- paste(setdiff(SET2, intersect(SET1,SET2)), collapse="\n")


setdiff1 <- t(as.data.frame(setdiff(SET1, intersect(SET1, SET2))))
setdiff2 <- t(as.data.frame(setdiff(SET2, intersect(SET1, SET2))))

#also no negative z score and filtered via LR table
write.csv(setdiff1, "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_significant_TF_PMF_cluster.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_significant_TF_PMF_cluster.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_significant_TF_PMF_cluster.pdf")
grid.draw(v)
dev.off()
##########################################################
#intracellular network
#ctrl
####
library(VennDiagram)
library(data.table)

csv1 <- fread("R_intra_network_ctrl.csv")
#csv1 <- results@intracellular_network_condition[["control"]]
csv2 <- data.table::fread("py_intra_network_ctrl.csv")

csv1 <- csv1[csv1$TF_Score > 0,]
csv2 <- csv2[csv2$TF_Score > 0,]
head(csv1)
head(csv2)
csv1 <- csv1[, !c("V1","TF_Score")]
#csv1$TF_Score <- NULL
csv2 <- csv2[, !c("V1","TF_Score")]


csv1_list <- apply(csv1, 1, paste, collapse = ",")
csv2_list <- apply(csv2, 1, paste, collapse = ",")
head(csv1_list, n=50)
head(csv2_list)

SET1 <- csv1_list 
SET2 <- csv2_list


v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY"),
  filename = NULL,
  output = TRUE
)

grid.newpage()
grid.draw(v)

setdiff1 <- setdiff(SET1, SET2)
setdiff2 <- setdiff(SET2, SET1)

fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_set_diff_intra_network_control_condition.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_set_diff_intra_network_control_condition.csv", col.names = FALSE)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_intra_network_diff_control_condition.pdf")
grid.draw(v)
dev.off()

####################################################################################
#intra network
#PMF MF2
########

csv1 <- fread("R_intra_network_PMF.csv")
#csv1 <- results@intracellular_network_condition[["PMF_MF2"]]
csv2 <- fread("py_intra_network_PMF.csv")

csv1 <- csv1[, !c("V1","TF_Score")]
#csv1$TF_Score <- NULL

csv2 <- csv2[, !c("V1","TF_Score")]

csv1_list <- apply(csv1, 1, paste, collapse = ",")
csv2_list <- apply(csv2, 1, paste, collapse = ",")

SET1 <- csv1_list 
SET2 <- csv2_list


v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY"),
  filename = NULL,
  output = TRUE
)


grid.newpage()
grid.draw(v)


setdiff1 <- setdiff(SET1, SET2)
setdiff2 <- setdiff(SET2, SET1)


fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_set_diff_intra_network_PMF_condition.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_set_diff_intra_network_PMF_condition.csv", col.names = FALSE)


pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_intra_network_diff_PMF_condition.pdf")
grid.draw(v)
dev.off()

##########################################################################
#intra network CLUSTER
#CONTROL

library(VennDiagram)
library(data.table)

csv1 <- fread("R_intra_network_ctrl_cluster.csv")
#csv1 <- results@intracellular_network_cluster[["control"]]
csv2 <- data.table::fread("py_intra_network_ctrl_cluster.csv")

csv1 <- csv1[csv1$TF_Score > 0,]
csv2 <- csv2[csv2$TF_Score > 0,]
head(csv1)
head(csv2)
csv1 <- csv1[, !c("V1","TF_Score")]
#csv1$TF_Score <- NULL
csv2 <- csv2[, !c("V1","TF_Score")]


csv1_list <- apply(csv1, 1, paste, collapse = ",")
csv2_list <- apply(csv2, 1, paste, collapse = ",")
head(csv1_list, n=50)
head(csv2_list)

SET1 <- csv1_list 
SET2 <- csv2_list


v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY"),
  filename = NULL,
  output = TRUE
)

grid.newpage()
grid.draw(v)

setdiff1 <- setdiff(SET1, SET2)
setdiff2 <- setdiff(SET2, SET1)

fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_set_diff_intra_network_ctrl_cluster.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_set_diff_intra_network_ctrl_cluster.csv", col.names = FALSE)

pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_intra_network_diff_ctrl_cluster.pdf")
grid.draw(v)
dev.off()

####################################################################################
#intra network CLUSTER
#PMF MF2
########

csv1 <- fread("R_intra_network_PMF_cluster.csv")
#csv1 <- results@intracellular_network_cluster[["PMF_MF2"]]
csv2 <- fread("py_intra_network_PMF_cluster.csv")

csv1 <- csv1[, !c("V1","TF_Score")]
#csv1$TF_Score <- NULL

csv2 <- csv2[, !c("V1","TF_Score")]

csv1_list <- apply(csv1, 1, paste, collapse = ",")
csv2_list <- apply(csv2, 1, paste, collapse = ",")

SET1 <- csv1_list 
SET2 <- csv2_list


v <- venn.diagram(
  x = list(SET1, SET2),
  category.names = c("Set R", "Set PY"),
  filename = NULL,
  output = TRUE
)


grid.newpage()
grid.draw(v)


setdiff1 <- setdiff(SET1, SET2)
setdiff2 <- setdiff(SET2, SET1)


fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/decoupler_main_scaled/R_set_diff_intra_network_PMF_cluster.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/decoupler_main_scaled/PY_set_diff_intra_network_PMF_cluster.csv", col.names = FALSE)


pdf("Venn_Diagrams_and_csvs/decoupler_main_scaled/venn_diagram_PY_R_intra_network_diff_PMF_cluster.pdf")
grid.draw(v)
dev.off()
