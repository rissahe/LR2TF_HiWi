#Main function that runs the IntraTalker Analysis
#Performs intracellular network anaysis based on conditions and clusters using decoupleR TF activity results and scRNA-seq data.

#Parameters:
#anndataobject: Input Anndata Object as h5ad file
#tf_activities: Matrix with TF activities for each cell in the scRNA-seq data; Input as either csv or DataFrame
#arguments_list: named list with custom options for the analysis

#import anndata as ad
#import pandas as pd
#import scanpy as sc
#import os
#import re 

def IntraTalker_analysis(anndataobject, tf_activities = None, arguments_list = None):
    
    #load Anndata object
    if (isinstance(anndataobject, str)):
        anndataobject = ad.read_h5ad(anndataobject)

    #get arguments list
    arguments_list = validate_input_arguments(arguments_list)

    #if the R version of decoupler was used, the tf activity matrix needs to be transposed
    if arguments_list["decoupler_matrix_format"] == "R":
        anndataobject = anndataobject.T

    #create directory
    if not os.path.isdir(arguments_list["out_path"]):
        os.mkdir(arguments_list["out_path"])
        tf_path = arguments_list["out_path"] + "TF_results/"
        os.mkdir(tf_path)
    else:
        tf_path = arguments_list["out_path"] + "TF_results/"


    condition = anndataobject.obs[arguments_list["condition"]]

    #load tf activities(csv input)
    if isinstance(tf_activities, str):
         tf_activities = ad.read_csv(tf_activities)
         tf_activities.obs = anndataobject.obs.reindex(tf_activities.obs.index)
         tf_activities.obsm = anndataobject.obsm
         tf_activities.uns = anndataobject.uns

    #load tf activities(DataFrame input)
    elif isinstance(tf_activities, pd.DataFrame):
         tf_activities = ad.AnnData(tf_activities)
         tf_activities.obs = anndataobject.obs.reindex(tf_activities.obs.index)
         tf_activities.obsm = anndataobject.obsm
         tf_activities.uns = anndataobject.uns

    elif tf_activities is None:
         raise NameError("Please attach a csv file with the tf activity values. (For further clarification view the 'Decoupler' section of the vignette.)")
    
    
    sc.pp.scale(tf_activities)

    #sets the stage for decision if single condition or comparison analysis is done
    if not arguments_list["comparison_list"] is None:
        if (len(arguments_list["comparison_list"]) > 0) & (len(pd.unique(anndataobject.obs[arguments_list["condition"]])) < 2):
            arguments_list["comparison_list"] = None
            print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")

    #code for only single condition analysis/cluster analysis where marker TFs are grouped by celltype/cluster
    if arguments_list["comparison_list"] is None:

        result_list = {}
        gene_expression_list = {}
        CTR_cluster_list = {}
        intranet_cluster_list = {}

        #loops through each condition and creates single condition results
        for name_iterable in anndataobject.obs[arguments_list["condition"]].unique():
            #subset anndata object and tf activities to single condition
            sub_object = anndataobject[anndataobject.obs[arguments_list["condition"]] == name_iterable].copy()
            tf_activities_sub = tf_activities[tf_activities.obs[arguments_list["condition"]] == name_iterable].copy()

            #get average gene expression for single condition from scRNA data
            sub_object_avg = AverageExpression(sub_object, name_iterable= name_iterable, celltype = arguments_list["celltype"], outpath= arguments_list["out_path"])
        
            #get marker TFs for single condition
            tf_activity_scores = get_significant_tfs(tf_activities_sub,
                                                        name_iterable,
                                                        tf_path,
                                                        None,
                                                        celltype = arguments_list["celltype"],
                                                        pval = arguments_list["pval"],
                                                        meanchange = arguments_list["meanchange"],
                                                        plot = arguments_list["plot"],
                                                        condition_comparison= False)
            
            print("tf activities done")

            result_list[name_iterable] = tf_activity_scores
            gene_expression_list[name_iterable] = sub_object_avg
            
    
            CTR_cluster_list[name_iterable] = generate_CrossTalkeR_input(tf_activity_scores["cluster"],
                                                                            gene_expression_list[name_iterable],
                                                                            arguments_list["out_path"],                                             
                                                                            arguments_list["reg"],
                                                                            arguments_list["organism"])

            print("CTR input cluster done")
              
            intranet_cluster_list[name_iterable] = generate_intracellular_network(tf_activity_scores["cluster"],
                                                                                  gene_expression_list[name_iterable],
                                                                                  arguments_list["out_path"],
                                                                                  arguments_list["reg"],
                                                                                  arguments_list["organism"])
            
            print("intracellular network cluster done")
            
        #save as result object
        tf = make_TFOBj(
                tf_activities_condition = list(),
                tf_activities_cluster = result_list,
                average_gene_expression = gene_expression_list,
                regulon = arguments_list["reg"],
                CTR_input_condition = list(),
                CTR_input_cluster = CTR_cluster_list,
                intracellular_network_condition = list(),
                intracellular_network_cluster = intranet_cluster_list)

        #with open((arguments_list["out_path"] + "tf.pickle"), "wb") as file:
        #    pickle.dump(tf, file)

        return tf

    else:
        #code for compared condition analysis with additional single condition analysis
        
        #compared condition analysis:
        out_path_compared = (tf_path + "compared")
        if not os.path.isdir(out_path_compared):
            os.mkdir(out_path_compared)

        #finds marker TFs for compared conditions
        compared_significant_tfs = condition_comparison_significant(tf_activities, out_path_compared, arguments_list["celltype"], 
                                                                    arguments_list["condition"], arguments_list["comparison_list"], 
                                                                    arguments_list["num_cell_filter"])
        
        print("compared tfs done")
        
        #creates heatmaps for compared condition (TF activity in first condition when compared to second condition)
        if arguments_list["plot"] == True:
            plot_condition_tf_activities(compared_significant_tfs, out_path_compared)
            plot_condition_tf_activities_compressed(compared_significant_tfs, out_path_compared)

    
        result_condition_list = {}
        result_cluster_list = {}
        gene_expression_list = {}
        CTR_condition_list = {}
        CTR_cluster_list = {}
        intranet_condition_list = {}
        intranet_cluster_list = {}

        #single condition/cluster analysis:
        for name_iterable in anndataobject.obs[arguments_list["condition"]].unique():
            #subset anndata object and tf activities to single condition
            sub_object = anndataobject[anndataobject.obs[arguments_list["condition"]] == name_iterable]
            tf_activities_sub = tf_activities[tf_activities.obs[arguments_list["condition"]] == name_iterable]
            
            #initialize df
            compared_tfs = pd.DataFrame({"gene" : pd.Series(dtype="str"), "tag" : pd.Series(dtype="str"), "cluster" : pd.Series(dtype="str")})
        
            #filters compared tfs
            for result_name, df in compared_significant_tfs.items(): 
                if name_iterable in result_name:
                    tf_condition_significant = compared_significant_tfs[result_name]
                    tf_condition_significant = tf_condition_significant[tf_condition_significant["FDR"] < arguments_list["pval"]]
                    tf_condition_significant = tf_condition_significant[(tf_condition_significant["meanchange"] > float(arguments_list["meanchange"])) | (tf_condition_significant["meanchange"] < (0 - float(arguments_list["meanchange"])))]
                    tf_condition_significant = tf_condition_significant[["tf", "tag", "CellType"]]
                    tf_condition_significant.rename(columns={"tf":"gene", "CellType": "cluster"}, inplace=True)
                    compared_tfs = pd.concat([compared_tfs, tf_condition_significant])

            #replaces special symbols in condition names with _
            re.sub("([,;.:-])", "_", name_iterable)

            #get average gene expression for single condition from scRNA data
            sub_object_avg = AverageExpression(sub_object, name_iterable= name_iterable, celltype = arguments_list["celltype"], 
                                               outpath= arguments_list["out_path"])
            
            #get marker TFs for single condition
            tf_activity_scores = get_significant_tfs(tf_activities_sub,
                                               name_iterable,
                                               tf_path,
                                               compared_tfs,
                                               celltype = arguments_list["celltype"],
                                               pval = arguments_list["pval"],
                                               meanchange = arguments_list["meanchange"],
                                               plot = arguments_list["plot"],
                                               condition_comparison= True)
            
            print("tf activities done")

            result_condition_list[name_iterable] = tf_activity_scores["condition"]
            result_cluster_list[name_iterable] = tf_activity_scores["cluster"]
            gene_expression_list[name_iterable] = sub_object_avg
            
        
            CTR_condition_list[name_iterable] = generate_CrossTalkeR_input(tf_activity_scores["condition"],
                                                                            gene_expression_list[name_iterable],
                                                                           arguments_list["out_path"],                                             
                                                                           arguments_list["reg"],
                                                                           arguments_list["organism"])
            
            print("CTR input condition done")
    
            CTR_cluster_list[name_iterable] = generate_CrossTalkeR_input(tf_activity_scores["cluster"],
                                                                            gene_expression_list[name_iterable],
                                                                            arguments_list["out_path"],                                             
                                                                            arguments_list["reg"],
                                                                            arguments_list["organism"])
    
            print("CTR input cluster done")

            intranet_condition_list[name_iterable] = generate_intracellular_network(tf_activity_scores["condition"],
                                                                                  gene_expression_list[name_iterable],
                                                                                  arguments_list["out_path"],
                                                                                  arguments_list["reg"],
                                                                                  arguments_list["organism"])
    
            print("ntracellular network condition done")

            intranet_cluster_list[name_iterable] = generate_intracellular_network(tf_activity_scores["cluster"],
                                                                                  gene_expression_list[name_iterable],
                                                                                  arguments_list["out_path"],
                                                                                  arguments_list["reg"],
                                                                                  arguments_list["organism"])
            print("intracellular network cluster done")
        
        #save as result object
        tf = make_TFOBj(
                tf_activities_condition = result_condition_list,
                tf_activities_cluster = result_cluster_list,
                average_gene_expression = gene_expression_list,
                regulon = arguments_list["reg"],
                CTR_input_condition = CTR_condition_list,
                CTR_input_cluster = CTR_cluster_list,
                intracellular_network_condition = intranet_condition_list,
                intracellular_network_cluster = intranet_cluster_list)

        #with open((arguments_list["out_path"] + "tf.pickle"), "wb") as file:
        #    pickle.dump(tf, file)
        
        return tf