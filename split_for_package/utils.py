#' Add entry to a dataframe in the CrossTalkeR input format
#'
#' @export
def add_entry(source, target, gene_A, gene_B, type_gene_A, type_gene_B, MeanLR):
  df = {"source" : source,
      "target" : target,
      "gene_A" : gene_A,
      "gene_B" : gene_B,
      "type_gene_A" : type_gene_A,
      "type_gene_B" : type_gene_B,
      "MeanLR" : MeanLR}

  return df

def map_t_value(tf_scores_df, anndataobject_markers):
    genes = []
    clusters = []
    t_values = []  
    anndataobject_markers = anndataobject_markers.set_index("gene")
    
    for i in range(len(anndataobject_markers.index)):
        a = anndataobject_markers.index[i]
        for j in range(len(tf_scores_df.index)):
            b = tf_scores_df.index[j]
            if a == b:
                c = anndataobject_markers.columns.get_loc("cluster")
                cluster = anndataobject_markers.iloc[i, c]
                gene_rows = tf_scores_df.iloc[j]
                if isinstance(gene_rows, pd.Series):
                    gene_rows = gene_rows.to_frame().T

                for _, gene_row in gene_rows.iterrows():
                     score = gene_row[cluster]
                     genes.append(a)
                     clusters.append(cluster)
                     t_values.append(score)

    t_value_df = pd.DataFrame({
        'gene': genes,
        'cluster': clusters,
        't_value': t_values
    })

    return t_value_df

def create_unfiltered_tf_scores(tf_scores_df, condition, celltype, out_path):   
    summarized_tf_scores_df = tf_scores_df.groupby(celltype, observed = True).mean().T
    #tf_scores_df.groupby(celltype, observed = True).apply(display)
    #agg(["mean", "var"])
    summarized_tf_scores_df.to_csv(out_path + "/unfiltered_tf_scores_" + condition + ".csv")
    return summarized_tf_scores_df

########################################################################################
#' Saves most variable transcription factor activity scores per cell type into a table
#'
#' This function saves the transcription factor activity scores per cell type
#' into a csv table.
#'
#' @param tf_scores data frame with transcription factor activity scores per cell type
#' @param condition Sample condition for file naming(e.g. control, disease ...)
#' @param out_path Output path to save results
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import stringr
#' @export
def save_variable_tf_score(filtered_summarized_tf_scores_df, condition, out_path, plot):
    filtered_summarized_tf_scores_df["var"] = filtered_summarized_tf_scores_df.apply(statistics.variance, axis=1).unique()
    filtered_summarized_tf_scores_df.to_csv(out_path + "/variable_tf_scores_" + condition + ".csv")

    if plot:
        top_variable_tfs = filtered_summarized_tf_scores_df.sort_values("var", ascending=False).head(n=20).drop(columns="var")
        cluster_map = sns.clustermap(top_variable_tfs, cmap="vlag", center=0, vmin=top_variable_tfs.min(axis=None), cbar_kws={"label": "t-score"},
                                    cbar_pos=(1, 0.5, 0.02, 0.1), dendrogram_ratio=(0.2, 0.05))
        cluster_map.ax_heatmap.set_xlabel("Cell Type")
        cluster_map.ax_heatmap.set_ylabel("Transcription Factor")
        plt.setp(cluster_map.ax_heatmap.get_xticklabels(), rotation=45) 
        plt.savefig(out_path + "/tf_activity_top20_variable_" + condition + ".pdf", bbox_inches='tight')
        plt.close()

    filtered_summarized_tf_scores_df_var = filtered_summarized_tf_scores_df
    return filtered_summarized_tf_scores_df_var

def eval_pval(p_val):
    p_val = float(p_val)
    if p_val < 0.001: 
      txt = "***"
    elif p_val < 0.01: 
      txt = "**"
    elif p_val < 0.05: 
      txt = "*"
    else:
      txt = "ns"
    return(txt)


def eval_meanchange_tag(meanchange):
    if meanchange >= 1.0: 
      txt = "***"
    elif meanchange > 0.5: 
      txt = "**"
    elif meanchange > 0.0: 
      txt = "*"
    else:
      txt = "ns"
    return(txt)

def AverageExpression(sub_object, celltype = None, name_iterable = None, outpath = None):
    gene_ids = sub_object.var.index.values
    obs = sub_object[:,gene_ids].X.toarray()
    obs = np.expm1(obs)
    avg_df = pd.DataFrame(obs,columns=gene_ids,index= sub_object.obs[celltype])
    avg_df = avg_df.groupby(level=0, observed=False).mean()
    avg_df.T.to_csv(outpath + name_iterable + "_average_gene_expression_by_cluster_exp.csv")

    return avg_df.T

#############################################################################################
#' Check arguments passed by User for validity
#'
#' @param arguments_list list of user defined arguments
#' @return val_arguments validated list of arguments
#' @export
def validate_input_arguments (arguments_list):
    if arguments_list["out_path"] is None:
        print("Please provide an output path")
    elif arguments_list["out_path"][-1] != "/":
        arguments_list["out_path"] = arguments_list["out_path"] + "/"

    if arguments_list["celltype"] is None:
        print("Please provide the name of the metadata field containing cell type annotations")

    if arguments_list["condition"] is None:
        print("Please provide the name of the metadata field containing condition annotations")

    if arguments_list["organism"] is None:
        arguments_list["organism"] = "human"

    #if arguments_list["comparison_list"] is None:
    #    arguments_list["comparison_list"] = np.nan

    if arguments_list["meanchange"] is None:
        arguments_list["meanchange"] = 0.0

    if arguments_list ["pval"] is None:
        arguments_list["pval"] = 0.05

    if arguments_list ["num_cell_filter"] is None:
        arguments_list["num_cell_filter"] = 0

    if arguments_list["reg"] is None:
        raise ValueError("Please provide a regulon csv.")

    elif isinstance(arguments_list["reg"], str):
        arguments_list["reg"] = pd.read_csv(arguments_list["reg"], index_col=0)
        arguments_list["reg"] = pd.DataFrame.rename(arguments_list["reg"], columns={"source" : "tf"})

    #is the naming of tf fine like this?
    if not "tf" in arguments_list["reg"] and "target" in arguments_list["reg"] and "weight" in arguments_list["reg"]:
        raise NameException("Not all necessary columns found in regulon table! Please make sure that the regulon has the columns 'source', 'target' and 'weight'!")
    
    if arguments_list["plot"] is None:
        arguments_list["plot"] = True
    elif not isinstance(arguments_list["plot"], (bool)):
        raise ValueError("lot argument must be a boolean value.")
        
   
    return(arguments_list)


##############################################################################################################
#' Adds gene type to the name of the gene
#'
#' Description
#'
#' @param df dataframe with all interactions
def add_node_type(df):
    
    df['gene_A'] = np.where(df['type_gene_A'] == 'Ligand', df['gene_A'] + '|L', df['gene_A'])
    df['gene_A'] = np.where(df['type_gene_A'] == 'Receptor', df['gene_A'] + '|R', df['gene_A'])
    df['gene_A'] = np.where(df['type_gene_A'] == 'Transcription Factor', df['gene_A'] + '|TF', df['gene_A'])
    df['gene_B'] = np.where(df['type_gene_B'] == 'Ligand', df['gene_B'] + '|L', df['gene_B'])
    df['gene_B'] = np.where(df['type_gene_B'] == 'Receptor', df['gene_B'] + '|R', df['gene_B'])
    df['gene_B'] = np.where(df['type_gene_B'] == 'Transcription Factor', df['gene_B'] + '|TF', df['gene_B'])
    return df


##############################################################################################################
#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_prediction path to or dataframe with ligand-receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
def combine_LR_and_TF(tf_table, LR_prediction, out_path, condition, add_nodetype = False):

  if isinstance(LR_prediction, pd.DataFrame):
    lr_table = LR_prediction
  else: 
    lr_table = pd.read_csv(LR_prediction, index_col=0)
  

  intra_connections = pd.DataFrame()
  for celltype in np.unique([lr_table["source"], lr_table["target"]]):
    lr_filtered_ligands = lr_table[lr_table["source"] == celltype]
    lr_filtered_receptors = lr_table[lr_table["target"] == celltype]
    lr_ligands = np.unique(lr_filtered_ligands["gene_A"])
    lr_receptors = np.unique(lr_filtered_receptors["gene_B"])
    #print(condition, celltype, lr_receptors)

    tf_table_receptors = tf_table[(tf_table["target"] == celltype) & (tf_table["type_gene_A"] == "Receptor")]
    tf_table_ligands = tf_table[(tf_table["source"] == celltype) & (tf_table["type_gene_B"] == "Ligand")]

    tf_receptor_interactions =  tf_table_receptors[tf_table_receptors["gene_A"].isin(lr_receptors)]
    tf_ligand_interactions = tf_table_ligands[tf_table_ligands["gene_B"].isin(lr_ligands)]


    intra_connections = pd.concat([intra_connections, tf_receptor_interactions, tf_ligand_interactions])
  intra_connections["all_pair"] = (intra_connections["source"] + "/" 
                                    + intra_connections["gene_A"] + "/"
                                    + intra_connections["target"] + "/"
                                    + intra_connections["gene_B"])
    
  intra_connections = intra_connections.drop_duplicates(subset=["all_pair"])
  intra_connections.drop(columns=["all_pair"], inplace=True)

  complete_interactions = pd.concat([intra_connections, lr_table])

  if add_nodetype:
    complete_interactions = add_node_type(complete_interactions)
      
  complete_interactions.to_csv((out_path + "CrossTalkeR_input_" + condition + ".csv"), index=0)
  return(complete_interactions)


#######################################################################################################
#' Combining Ligand-Receptor interaction prediction with Transcription Factor interaction predictions considering receptor complexes
#'
#' Description
#'
#' @param tf_table table with tf interactions
#' @param LR_prediction path to or dataframe with ligand-receptor interaction prediction
#' @param out_path path to save results
#' @param condition sample condition of data
#' @return complete_interactions table with all interactions
def combine_LR_and_TF_complexes(tf_table, LR_prediction, out_path, condition, add_nodetype = False):

  if isinstance(LR_prediction, pd.DataFrame):
    lr_table = LR_prediction
  else: 
    lr_table = pd.read_csv(LR_prediction, index_col=False)
  
  intra_connections = pd.DataFrame()
  for celltype in np.unique([lr_table["source"], lr_table["target"]]):
    lr_filtered_ligands = lr_table[lr_table["source"] == celltype]
    lr_filtered_receptors = lr_table[lr_table["target"] == celltype]
    lr_ligands = np.unique(lr_filtered_ligands["gene_A"])
    lr_receptors = np.unique(lr_filtered_receptors["gene_B"])
    #print(lr_filtered_receptors)

    lr_receptors = pd.Series(lr_receptors)
    contains_complex = lr_receptors.str.contains("_", na=False)
    
    R_with_complex = lr_receptors[contains_complex]
    #print("R_with_complex", R_with_complex)
    R_without_complex = lr_receptors[(~contains_complex)]
  
    tf_table_receptors = tf_table[(tf_table["target"] == celltype) & (tf_table["type_gene_A"] == "Receptor")]
    tf_receptor_interactions =  tf_table_receptors[tf_table_receptors["gene_A"].isin(R_without_complex)]

    c_receptors = tf_table_receptors[tf_table_receptors["gene_A"].apply(lambda x: any(gene in x.split("+") for gene in lr_receptors))]
    #print("c receptors", c_receptors)

    complex_df = pd.DataFrame()
    if len(R_with_complex) > 0:
      for complex in R_with_complex:
        receptors = complex.split("_")
        R_TF_with_complex = tf_table_receptors[tf_table_receptors["gene_A"].isin(receptors)]
 
        if len(R_TF_with_complex) == 0:
          continue
        
        R_TF_with_complex.drop_duplicates()
        R_TF_with_complex.loc[:,"gene_A"] = complex
        complex_df = pd.concat([complex_df, R_TF_with_complex])

      #complex_df.drop_duplicates()

    tf_receptor_interactions = pd.concat([tf_receptor_interactions, complex_df])
    #print("tf_receptor_interactions", tf_receptor_interactions)
    
    tf_table_ligands = tf_table[(tf_table["source"] == celltype) & (tf_table["type_gene_B"] == "Ligand")]

    
    tf_ligand_interactions = tf_table_ligands[tf_table_ligands["gene_B"].isin(lr_ligands)]

    intra_connections = pd.concat([intra_connections, tf_receptor_interactions, c_receptors, tf_ligand_interactions])
  
  intra_connections["all_pair"] = (intra_connections["source"] + "/" 
                                    + intra_connections["gene_A"] + "/"
                                    + intra_connections["target"] + "/"
                                    + intra_connections["gene_B"])
  #print("intra_connections", intra_connections)
  intra_connections = intra_connections.drop_duplicates(subset=["all_pair"])
  intra_connections.drop(columns=["all_pair"], inplace=True)

  complete_interactions = pd.concat([intra_connections, lr_table])
  
  if add_nodetype:
    complete_interactions = add_node_type(complete_interactions)
      
  complete_interactions.to_csv((out_path + "CrossTalkeR_input_" + condition + ".csv"), index=False, quoting= csv.QUOTE_NONNUMERIC)
  return(complete_interactions)
