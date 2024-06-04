
import scanpy as sc
import decoupler as dc
import pandas as pd
import anndata as ad
import pooch

ann_data = sc.read_h5ad("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/anndata_object.h5ad")
reg = pd.read_csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/filterd_regulon.csv")

dc.run_ulm( mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

estimates =ann_data.obsm['ulm_estimate']
estimates.to_csv("/home/larissa/Documents/Larissa_HiWi/LR2TF_test/decoupler_results.csv")