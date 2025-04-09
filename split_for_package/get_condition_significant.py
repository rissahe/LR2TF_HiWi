#decoupler condition comparison
def condition_comparison_significant(tf_activities, out_path, celltype, condition, comparison_list, num_cell_filter = 0):

    vs_df_dic = {}

    if isinstance(comparison_list[0], str):  
        comparison_list = [comparison_list]
    
    for vs1, vs2 in comparison_list:

        print(f"vs1: {vs1}, vs2: {vs2}") 

        all_tf_list = tf_activities.var_names

        res = pd.DataFrame()
        for i in tf_activities.obs[celltype].unique(): 
            comparison_sub = tf_activities[(tf_activities.obs[celltype] == i) & (tf_activities.obs[condition].isin([vs1, vs2]))]
            if len(pd.unique(comparison_sub.obs[condition])) == 2:
                condition_table = comparison_sub.obs[[condition]].copy()
                condition_table.columns = ["condition"]
                metadata_counts = condition_table.groupby("condition", observed = False).size()
                
                if (metadata_counts.iloc[0] + metadata_counts.iloc[1]) > num_cell_filter:
                    g = comparison_sub.obs[condition].astype("category")
                    g = g.cat.set_categories([vs1, vs2])
        
                    #finding markers for the condition comparison
                    res_tmp = dc.rank_sources_groups(comparison_sub, groupby= condition, reference="rest", method="t-test")
                    res_tmp.rename(columns={"group": "condition", "statistic" : "scores"}, inplace=True)
            
                    #calculating wilcoxon scores to use for r value calculation (used in heatmaps)
                    res_heatmap = dc.rank_sources_groups(comparison_sub, groupby= condition, reference="rest", method="wilcoxon")
                    res_heatmap.rename(columns={"group": "condition", "statistic" : "scores"}, inplace=True)

                    group1 = comparison_sub.X[g == vs1]
                    group2 = comparison_sub.X[g == vs2]
                        
                    res_heatmap["r"] = (res_heatmap["scores"] / np.sqrt(len(group1) + len(group2)))
                    res_heatmap["CellType"] = i
                    res_tmp["CellType"] = i
                    _, res_tmp["FDR"], _, _ = multipletests(res_tmp["pvals"], alpha=0.05, method='fdr_bh')
                    
                    
                    res_tmp = res_tmp[res_tmp["condition"] == vs1]

                    res_heatmap = res_heatmap[["names", "r", "CellType", "condition"]]
                    res_tmp = pd.merge(res_tmp, res_heatmap, on = ["names","CellType", "condition"])

                    res = pd.concat([res, res_tmp], ignore_index=True)

        res_df = res.dropna()

        def assign_significance_tag(fdr):
            if fdr < 0.001:
                return "***"
            elif fdr < 0.01:
                return "**"
            elif fdr < 0.05:
                return "*"
            else:
                return "ns"

        res_df = res_df.assign(tag=res_df["FDR"].apply(assign_significance_tag))

        res_df.rename(columns={"names":"tf", "group": "condition"}, inplace=True)
        res_df.to_csv(f"{out_path}/all_tfs_{vs1}_vs_{vs2}.csv", index=False)

        result_name = f"{vs1}_{vs2}"
        vs_df_dic[result_name] = res_df
    return vs_df_dic