#install.packages('VennDiagram')
library(VennDiagram)

#########
#CTRL
#########

csv1 <- read.csv("R_CRT_no_negative_score.csv")
row.names(csv1) <- NULL
csv1 <- csv1[c(1:587),]

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
csv2 <- csv2[c(1:314),]
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
write.csv(setdiff1, "Venn_Diagrams_and_csvs/R_unique_CTRL_decoupler_scaled.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/Py_unique_CTRL_decoupler_scaled.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

#filtered_with_LR_table_no_neg_score_
pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_CTR_input_CTRL_decoupler_scaled.pdf")
grid.draw(v)
dev.off()

################

#PMF
####################

csv1 <- read.csv("R_CRT_no_negative_score_PMF.csv")
row.names(csv1) <- NULL
csv1 <- csv1[c(1:623),]

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
csv2 <- csv2[c(1:327),]
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

write.csv(setdiff1, "Venn_Diagrams_and_csvs/R_unique_PMF_decoupler_scaled.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/Py_unique_PMF_decoupler_scaled.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_CTR_input_PMF_decoupler_scaled.pdf")
grid.draw(v)
dev.off()

########################################################
#type gene = TF
#CTRL
################

csv1 <- read.csv("LR2TF_test_run/results/TF_results/control/significant_condition_tf_results_control.csv")
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
write.csv(setdiff1, "Venn_Diagrams_and_csvs/R_set_diff_only_TF_ctrl_decoupler_scaled.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/PY_set_diff_only_TF_ctrl_decoupler_scaled.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_significant_TF_condition_ctrl_decoupler_scaled.pdf")
grid.draw(v)
dev.off()
########################################################
#type =TF
#PMF.MF2
#####################

csv1 <- read.csv("LR2TF_test_run/results/TF_results/PMF_MF2/significant_condition_tf_results_PMF_MF2.csv")
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
write.csv(setdiff1, "Venn_Diagrams_and_csvs/R_significant_TF_condition_PMF_scipy_decoupler_scaled.csv", row.names = FALSE)
write.csv(setdiff2, "Venn_Diagrams_and_csvs/PY_significant_TF_condition_PMF_scipy_decoupler_scaled.csv", row.names = FALSE)


grid.newpage()
grid.draw(v)

pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_set_diff_only_TF_PMF_decoupler_scaled.pdf")
grid.draw(v)
dev.off()

##########################################################
#intracellular network
#ctrl
####
library(VennDiagram)
library(data.table) 

csv1 <- fread("R_intra_network_ctrl.csv")
csv2 <- fread("py_intra_network_ctrl.csv")
head(csv1)
head(csv2)
csv1 <- csv1[, !c("V1","TF_Score")]
csv2 <- csv2[, !c("V1","TF_Score")]

csv1_list <- apply(csv1, 1, paste, collapse = ",")
csv2_list <- apply(csv2, 1, paste, collapse = ",")
head(csv1_list, n=500)
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

fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/R_set_diff_intra_network_ctrl_decoupler_scaled.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/PY_set_diff_intra_network_ctrl_decoupler_scaled.csv", col.names = FALSE)

pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_intra_network_diff_ctrl_decoupler_scaled.pdf")
grid.draw(v)
dev.off()

####################################################################################
#intra network
#PMF MF2
########

csv1 <- fread("R_intra_network_PMF.csv")
csv2 <- fread("py_intra_network_PMF.csv")

csv1 <- csv1[, !c("V1","TF_Score")]
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


fwrite(data.table(setdiff1), "Venn_Diagrams_and_csvs/R_set_diff_intra_network_PMF_decoupler_scaled.csv", col.names = FALSE)
fwrite(data.table(setdiff2), "Venn_Diagrams_and_csvs/PY_set_diff_intra_network_PMF_decoupler_scaled.csv", col.names = FALSE)


pdf("Venn_Diagrams_and_csvs/venn_diagram_PY_R_intra_network_diff_PMF_decoupler_scaled.pdf")
grid.draw(v)
dev.off()