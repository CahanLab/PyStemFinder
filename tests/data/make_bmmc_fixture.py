"""Pack the R stemFinder reference export into tests/data/bmmc_r_reference.h5ad.

Run export_bmmc_reference.R first, then:

    python make_bmmc_fixture.py <R export dir> <bmmc_competitor_results.csv> <out .h5ad>

The competitor scores are https://cnobjects.s3.amazonaws.com/stemFinder/bmmc_competitor_results.csv.

Contents (3427 Tabula Muris bone marrow cells x 91 mouse S/G2M cell cycle genes):
    X                 Seurat scale.data (float32; stemFinder only uses its sign)
    layers['data']    Seurat log-normalized data (float64, sparse)
    obsp['distances'] Seurat RNA_nn kNN graph (k = 59, each cell includes itself)
    obs               Phenotype, Ground_truth, R scores 'R_<method>_stemFinder[_raw]',
                      and the competitor 'ccat_invert' / 'CytoTRACE_invert' scores
    uns               'R_markers': the markers R used, R performance metrics, and R session versions
"""
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.sparse

export, competitor_csv, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]

obs = pd.read_csv(export / "obs.csv", index_col="cell")
obs.index.name = None
obs["Phenotype"] = obs["Phenotype"].astype("category")
obs = obs.rename(columns={c: f"R_{c}" for c in obs.columns if "stemFinder" in c})
competitor = pd.read_csv(competitor_csv, index_col=0).loc[obs.index]
obs[["ccat_invert", "CytoTRACE_invert"]] = competitor[["ccat_invert", "CytoTRACE_invert"]]

markers = (export / "markers.txt").read_text().split()
genes = list(dict.fromkeys(markers))
scale_data = pd.read_csv(export / "scale_data.csv", index_col=0).loc[obs.index, genes]
data = pd.read_csv(export / "data.csv", index_col=0).loc[obs.index, genes]

knn = pd.read_csv(export / "knn.csv")
n = len(obs)
graph = scipy.sparse.csr_matrix((np.ones(len(knn), dtype=np.float32), (knn["i"], knn["j"])), shape=(n, n))

performance = pd.read_csv(export / "performance.csv").set_index("metric")["value"]
session = dict(line.split(" ", 1) for line in (export / "session.txt").read_text().splitlines())

adata = anndata.AnnData(
    X=scale_data.to_numpy(dtype=np.float32),
    obs=obs,
    var=pd.DataFrame(index=genes),
    layers={"data": scipy.sparse.csr_matrix(data.to_numpy(dtype=np.float64))},
    obsp={"distances": graph},
    uns={"R_markers": markers, "R_performance": performance.to_dict(), "R_session": session},
)
adata.write_h5ad(out, compression="gzip")
print(adata)
