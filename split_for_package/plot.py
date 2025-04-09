def h_clust(data):
            dist_matrix = pdist(data.T)
            linkage_matrix = linkage(dist_matrix, method = "average")   
            return linkage_matrix

def plot_condition_tf_activities(tf_activity_tables, out_path):
       
   for result_name, df in tf_activity_tables.items():
        name_df = pd.DataFrame(df)
        #print(name_df)
        significant_res = name_df[name_df["tag"] == "***"]
        significant_genes = np.unique(significant_res["tf"])
        name_df = name_df[name_df['tf'].isin(significant_genes)]
        #print(name_df)cha

        tag_mapping = name_df[["tf", "tag", "CellType"]]

        tag_mapping = tag_mapping.pivot(index="tf", columns="CellType", values="tag")
        tag_mapping.fillna("ns", inplace=True)

        cols_all_ns = tag_mapping.columns[(tag_mapping == "ns").all()]
        tag_mapping = tag_mapping.drop(columns=cols_all_ns)

    
        name_df_r = name_df[["r", "tf", "CellType"]]
        name_df_cluster = name_df_r.pivot_table(index = "tf", columns = "CellType", values = "r", aggfunc = "mean")
        name_df_cluster.fillna(0, inplace=True) 
        name_df_cluster = name_df_cluster.drop(columns=cols_all_ns)
        h_clust_matrix = h_clust(name_df_cluster)


      
        calc_size = ((len(name_df["CellType"].unique()) * 1.5), (len(name_df_cluster) * 0.2))
      
        cluster_map = sns.clustermap(name_df_cluster, cbar_kws={"label": "r"}, figsize=calc_size, cmap="vlag", center=0, annot= tag_mapping,
                      col_linkage= h_clust_matrix, fmt="", yticklabels=True, cbar_pos=(1, 0.5, 0.03, 0.05), dendrogram_ratio=(0.2, 0.05))
        plt.setp(cluster_map.ax_heatmap.get_xticklabels(), rotation=45)  

        plt.savefig(out_path + "/" + result_name + "_cluster_condition_activity_difference.pdf", bbox_inches='tight')
        plt.close()

def plot_condition_tf_activities_compressed(tf_activity_tables, out_path):
 
  
    for result_name, df in tf_activity_tables.items():
        name_df = pd.DataFrame(df)
        
        significant_res = name_df[name_df["tag"] == "***"]
        significant_genes = np.unique(significant_res["tf"])
        name_df = name_df[name_df['tf'].isin(significant_genes)]
        
        name_df_r = name_df[["r", "tf", "CellType"]]
        name_df_cluster = name_df_r.pivot_table(index = "tf", columns = "CellType", values = "r", aggfunc = "mean")
        name_df_cluster.fillna(0, inplace=True) 

        h_clust_matrix = h_clust(name_df_cluster)
        name_df_cluster.reset_index()

        calc_size = ((len(name_df["CellType"].unique()) * 1.5
        ), (len(name_df_cluster) * 0.06))
        cluster_map = sns.clustermap(name_df_cluster, cbar_kws={"label": "r"}, figsize=calc_size, cmap="vlag", center=0, annot= None,
                       yticklabels=False, col_linkage= h_clust_matrix, fmt="", cbar_pos=(1, 0.5, 0.02, 0.1), dendrogram_ratio=(0.2, 0.05))
        plt.setp(cluster_map.ax_heatmap.get_xticklabels(), rotation=45) 
        #cluster_map.figure.suptitle(result_name)
        
        
        plt.savefig(out_path + "/" + result_name +  "_cluster_condition_activity_difference_compressed.pdf", bbox_inches='tight')
        plt.close()


def plot_tf_activity(filtered_summarized_tf_scores_df, tag_mapping, condition, out_path):
    filtered_summarized_tf_scores_df = filtered_summarized_tf_scores_df.drop(columns="var")
    calc_size = ((len(filtered_summarized_tf_scores_df.columns.unique()) * 1.5), (len(filtered_summarized_tf_scores_df) * 0.06))
    cluster_map = sns.clustermap(filtered_summarized_tf_scores_df, cbar_kws={"label": "t-score"}, figsize=calc_size, cmap="vlag", center=0, 
                                yticklabels=False, cbar_pos=(1, 0.5, 0.02, 0.1), dendrogram_ratio=(0.2, 0.05))
    cluster_map.ax_heatmap.set_xlabel("Cell Type")
    cluster_map.ax_heatmap.set_ylabel("Transcription Factor")
    plt.setp(cluster_map.ax_heatmap.get_xticklabels(), rotation=45) 
    
    plt.savefig(out_path + "/tf_activity_compressed_" + condition + ".pdf", bbox_inches='tight')
    plt.close()


    cols_all_ns = tag_mapping.columns[(tag_mapping == "ns").all()]
    tag_mapping = tag_mapping.drop(columns=cols_all_ns)
    filtered_summarized_tf_scores_df = filtered_summarized_tf_scores_df.drop(columns=cols_all_ns)

    calc_size = ((len(filtered_summarized_tf_scores_df.columns.unique()) * 1.5), (len(filtered_summarized_tf_scores_df) * 0.2))
    cluster_map = sns.clustermap(filtered_summarized_tf_scores_df, cbar_kws={"label": "t-score"}, figsize=calc_size, cmap="vlag", center=0, annot= tag_mapping, fmt="", 
                                yticklabels=True, cbar_pos=(1, 0.5, 0.03, 0.05), dendrogram_ratio=(0.2, 0.05))
    cluster_map.ax_heatmap.set_xlabel("Cell Type")
    cluster_map.ax_heatmap.set_ylabel("Transcription Factor")
    plt.setp(cluster_map.ax_heatmap.get_xticklabels(), rotation=45) 

    plt.savefig(out_path + "/tf_activity_" + condition + ".pdf", bbox_inches='tight')
    plt.close()