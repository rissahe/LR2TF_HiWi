#' Analyse transcription factor activities for significant transcription factors
#'
#' Description
#'
#' @param tf_activities_sub: Input Anndata Object with TF activities as X matrix
#' @param condition: Experiment condition (e.g. disease, knockout ...)
#' @param out_path Output path to save results
#' @param tf_condition_significant: condition comparison results
#' @param celltype: variable that accesses celltype meta data in anndata object
#' @param pval p-value to filter results
#' @param meanchange: meanchange value to filter results
#' @param plot: boolean variable to decide if code plots heatmaps or not
#' @param condition_comparison: comparision between multiple conditions or single condition
#' @return A data frame with transcription factor activity scores per cell type

def get_significant_tfs(tf_activities_sub, condition, out_path, tf_condition_significant, celltype, pval, meanchange, plot, condition_comparison = False):
    
    renamed_condition = condition.replace(",", "_")

    #create directory
    single_result_path = out_path + renamed_condition 
    if not os.path.isdir(single_result_path):
        os.mkdir(single_result_path)

    number_of_clusters = len(tf_activities_sub.obs[celltype].cat.categories) 

    #get marker TFs for single condition 
    anndataobject_markers = dc.rank_sources_groups(tf_activities_sub, groupby= celltype, reference="rest", method="wilcoxon")
    anndataobject_markers.rename(columns={"names" : "gene", "group": "cluster", "statistic" : "scores"}, inplace=True)
    
    anndataobject_markers["tag"] = None
    anndataobject_markers["meanchange_tag"] = None
    
    #add significance tags to p-value and meanchange
    anndataobject_markers["tag"] = anndataobject_markers["pvals_adj"].apply(eval_pval)
    anndataobject_markers["meanchange_tag"] = anndataobject_markers["meanchange"].apply(eval_meanchange_tag)
     
    anndataobject_markers.to_csv(single_result_path + "/" + renamed_condition + "_specific_markers_t_test.csv",index=0)

    #using the custom pval and meanchange filters to create tag mapping for significant tfs
    clusters = sorted(anndataobject_markers["cluster"].unique())
    tag_mapping = anndataobject_markers[["gene", "tag", "meanchange_tag", "cluster", "pvals_adj", "meanchange"]]
    tag_mapping = tag_mapping[(tag_mapping["pvals_adj"] < float(pval))] 
    tag_mapping = tag_mapping[(tag_mapping["meanchange"] > float(meanchange)) | 
                              (tag_mapping["meanchange"] < -float(meanchange))]

    tag_mapping = tag_mapping.pivot(index="gene", columns="cluster", values="tag")

    for cluster in clusters:
        if cluster not in tag_mapping.columns:
            tag_mapping[cluster] = np.nan

    tag_mapping = tag_mapping[clusters]
    tag_mapping = tag_mapping.astype("object")
    tag_mapping.fillna("ns", inplace=True)

    #creating df with genes as index and tf activity scores summarized by celltype/cluster (columns)
    tf_activities_sub.obs_names = tf_activities_sub.obs[celltype].astype(str)
    tf_scores_df = tf_activities_sub.to_df()
    #tf_scores_df.columns.name = None
    unfiltered_tf_scores = create_unfiltered_tf_scores(tf_scores_df, condition, celltype, single_result_path)
   
    #Filter to only include tfs that match the tag_mapping/are markers
    col_num = tf_scores_df.columns.isin(tag_mapping.index)  
    filtered_tf_scores_df = tf_scores_df.loc[:, col_num]
    filtered_summarized_tf_scores_df = filtered_tf_scores_df.groupby(celltype, observed = False).mean().T
    filtered_summarized_tf_scores_df.sort_index(axis=1, inplace=True)
    filtered_summarized_tf_scores_df.index.name = "gene"
    filtered_summarized_tf_scores_df.to_csv(f"{single_result_path}/tf_scores_{condition}.csv")

    #includes variability value of a gene's tf activity over all celltypes/clusters
    tf_scores_variable_table = save_variable_tf_score(filtered_summarized_tf_scores_df, condition, single_result_path, plot)

    #plots heatmap
    if plot:
        plot_tf_activity(filtered_summarized_tf_scores_df, tag_mapping, condition, single_result_path)
    
    filtered_summarized_tf_scores_df.index = filtered_summarized_tf_scores_df.index.map(lambda x: re.sub(".,", "_", x))

    #merges the tf score df with the marker tf df so that only the tf scores of tfs classified as significant by decoupler rank sources groups function are included
    anndataobject_markers_merged = anndataobject_markers.merge(map_t_value(filtered_summarized_tf_scores_df, anndataobject_markers),  on=['gene', 'cluster'], how='inner')
    anndataobject_markers_merged = anndataobject_markers_merged[anndataobject_markers_merged.tag != "ns"]
    anndataobject_markers_merged = anndataobject_markers_merged.where(anndataobject_markers_merged["t_value"] > 0, np.nan) 
    anndataobject_markers_merged.dropna(inplace=True)
    res_tmp = anndataobject_markers_merged[["gene","tag", "cluster", "t_value"]]

    res_tmp.to_csv(single_result_path + "/significant_cluster_tf_results_" + renamed_condition + ".csv", index=0)
    res = {}
    res["cluster"] = res_tmp
    

    if condition_comparison:
        unfiltered_tf_scores = unfiltered_tf_scores.where(unfiltered_tf_scores > 0, np.nan) 
        tf_condition_significant["gene"] = tf_condition_significant["gene"].apply(lambda x: re.sub(".,", "_", x))
        tf_condition_significant = tf_condition_significant.merge(map_t_value(unfiltered_tf_scores, tf_condition_significant), left_on=None, right_on=None, left_index=False, right_index=False)
        tf_condition_significant.dropna(inplace=True)
    
        res["condition"] = tf_condition_significant
        res["condition"].to_csv(f"{single_result_path}/significant_condition_tf_results_{condition}.csv", index=0)

    return res