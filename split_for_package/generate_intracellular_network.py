def generate_CrossTalkeR_input(tf_activities, gene_expression, out_path, regulon = None, organism = "human"):

  if organism == "human":
    ligands = pd.read_csv("ligands_human.csv")
    R2TF = pd.read_csv("rtf_db_human.csv")
  elif organism == "mouse": 
    ligands = pd.read_csv("ligands_mouse.csv")
    R2TF = pd.read_csv("rtf_db_mouse.csv")
  else:
    NameError("Invalid organism to generate CrossTalkeR input!")

  ligands = ligands.drop_duplicates()
  R2TF = R2TF.set_index("tf")

  sorted_regulon = regulon[["tf", "target"]]
  sorted_regulon = sorted_regulon.set_index("tf")

  tf_activities = tf_activities[tf_activities["t_value"] > 0]
  output_list = []
  df_list_l = {}
  df_list_r = {}
  tf_activities.reset_index(drop = True)

  for row in range(len((tf_activities))):
    
    targets = []
    tf_ligands = []
    tf = str(tf_activities["gene"].iloc[row])

    if tf in sorted_regulon.index:
      targets = sorted_regulon.loc[tf, "target"]
      if isinstance(targets, str):
        targets = [targets]
      
    if tf in R2TF.index:
        receptors = R2TF.loc[tf, "receptor"]
        if isinstance(receptors, str):
          receptors = [receptors]
          
    else:
        receptors = []


    if len(targets) > 0:
      tf_ligands = np.intersect1d(targets, ligands)

    if organism == "human":
      if len(tf_ligands) > 0:
        existing_entries = set()
        for ligand in tf_ligands:
            expressed = False
            if ligand in gene_expression.index:
              ex_value = gene_expression.loc[ligand, tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")]]
              if (ex_value != 0):
                expressed = True

            if (expressed == True):
              df_list_l = add_entry(source = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                      target = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                      gene_A =tf_activities.iloc[row, tf_activities.columns.get_loc("gene")],
                                                      gene_B = ligand,
                                                      type_gene_A = "Transcription Factor",
                                                      type_gene_B= "Ligand",
                                                      MeanLR= tf_activities.iloc[row, tf_activities.columns.get_loc("t_value")]
                                                      )
              
              if (ligand, tf) not in existing_entries:
                existing_entries.add((ligand, tf))
                output_list.append(df_list_l)
              

      if (len(receptors) > 0):
        for receptor in receptors:
          df_list_r = add_entry(source =  tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                target = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                gene_A= receptor,
                                                gene_B= tf_activities.iloc[row, tf_activities.columns.get_loc("gene")],
                                                type_gene_A= "Receptor",
                                                type_gene_B= "Transcription Factor",
                                                MeanLR= tf_activities.iloc[row, tf_activities.columns.get_loc("t_value")]
                                                )
          output_list.append(df_list_r)

      
          
    else: 
      if len(tf_ligands) > 0:
        for ligand in tf_ligands:
            expressed = False
            for l in ligands:
              if l in gene_expression.index:
                ex_value = gene_expression.loc[l, tf_activities.iloc[row, 2]]
                if (ex_value != 0):
                  expressed = True
          
              df_list_l = {}
            
              if (expressed == True):
                df_list_l = add_entry(source = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                      target = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                      gene_A =tf_activities.iloc[row, tf_activities.columns.get_loc("gene")],
                                                      gene_B = ligand,
                                                      type_gene_A = "Transcription Factor",
                                                      type_gene_B= "Ligand",
                                                      MeanLR= tf_activities.iloc[row, tf_activities.columns.get_loc("t_value")]
                                                      )
                
                output_list.append(df_list_l)
      
      if (len(receptors) > 0):
        for receptor in receptors:
          df_list_r = {}
          df_list_r = add_entry(source =  tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                target = tf_activities.iloc[row, tf_activities.columns.get_loc("cluster")],
                                                gene_A= receptor,
                                                gene_B= tf_activities.iloc[row, tf_activities.columns.get_loc("gene")],
                                                type_gene_A= "Receptor",
                                                type_gene_B= "Transcription Factor",
                                                MeanLR= tf_activities.iloc[row, tf_activities.columns.get_loc("t_value")]
                                                )
                                                
          output_list.append(df_list_r)
        #tf_l.to_csv(single_result_path + "/" + renamed_condition + "_ctr_test_l.csv", index=0)
        #r_tf.to_csv(single_result_path + "/" + renamed_condition + "_ctr_test_r.csv", index=0)
   

  output_df = pd.DataFrame(output_list)
  output_df["gene_A"] = output_df["gene_A"].apply(lambda x: re.sub("_", "+", x))
  output_df["gene_B"] = output_df["gene_B"].apply(lambda x: re.sub("_", "+", x))
  output_df.drop_duplicates(inplace=True)
  output_df.to_csv("tf_l_r_R_data_cond_contr.csv", index=0)


  return output_df



def generate_intracellular_network(tf_activities, gene_expression, outpath, regulon, organism="human"):

    if len(tf_activities.shape) > 0:
        if organism == "human":
            R2TF = pd.read_csv("rtf_db_human.csv").set_index("tf")
        else:
            R2TF = pd.read_csv("rtf_db_mouse.csv").set_index("tf")

    sorted_regulon = regulon[["tf", "target"]].set_index("tf")

    #preextract values
    tf_genes = tf_activities["gene"].values
    tf_celltypes = tf_activities.iloc[:, 2].values
    tf_scores = tf_activities.iloc[:, 3].values

    TFTG_list = []
    RTF_list = []

    tf_activities = tf_activities[tf_activities["t_value"] > 0]
    tf_activities.reset_index(drop = True)
    for row in range(len(tf_activities)):
        tf = str(tf_genes[row])
        celltype = tf_celltypes[row]
        tf_score = tf_scores[row]

        targets = sorted_regulon.loc[tf, "target"] if tf in sorted_regulon.index else []
        #print(targets)
        if tf in R2TF.index:
           receptors = R2TF.loc[tf, "receptor"]
           receptors = [receptors] if isinstance(receptors, str) else receptors
        else:
           receptors = []

        if len(targets) > 0 and len(receptors) > 0:
            for target in targets:
          
                if target in gene_expression.index:
                    ex_value = gene_expression.at[target, celltype]
         
                    if ex_value != 0:
                        TFTG_list.append({
                            "celltype": celltype,
                            "TF": tf,
                            "Target_Gene": target,
                            "TF_Score": tf_score
                        })

            for receptor in receptors:
                RTF_list.append({
                    "TF": tf,
                    "Receptor": receptor
                })

    TFTG_df = pd.DataFrame(TFTG_list)
    RTF_df = pd.DataFrame(RTF_list)

    recept_regulon = pd.merge(RTF_df, TFTG_df, on="TF")
    recept_regulon = recept_regulon.sort_values("TF")
    return recept_regulon