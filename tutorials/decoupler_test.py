
import scanpy as sc
import decoupler as dc
import pandas as pd

#import subprocess
#import sys
#subprocess.check_call([sys.executable, "-m", "pip", "install", "git+https://github.com/saezlab/decoupler-py"])

#anndata_object.h5ad for decoupler

ann_data = sc.read_h5ad("D:\\studium\\HiWi\\LR2TF\\anndata_object.h5ad")
reg = pd.read_csv("D:\\studium\\HiWi\\LR2TF\\filterd_regulon.csv")

dc.run_ulm( mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

estimates =ann_data.obsm['ulm_estimate']
estimates.to_csv("D:\\studium\\HiWi\\LR2TF\\decoupler_results.csv")