loom = "prefix.loom"
rds2loomfile = "rds2loomfile.txt" # 单细胞与loom名对照表
coord2d = "umap.txt" #单细胞UMAP坐标
Celllist = "Clusters.txt" # 细胞名对cluster表

min_shared_counts1 = 20 # 在剪接/未剪接矩阵共享的counts数最小值
n_top_genes1 = 2000 # 需要计算的Hvg个数
n_pcs1 = 30  # 使用的PCs
n_neighbors1 = 30 # findneighbor步骤使用的nneighbors

import loompy
import scvelo as scv
import scanpy as sc
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scvelo as scv
import itertools
import os
import sys
from pandas.api.types import CategoricalDtype


with loompy.connect(loom) as ds:
  CellID = list(ds.ca.CellID)
  Genename=list(ds.ra.Gene)

rds2loom = pd.read_csv(rds2loomfile,header= 0,sep="\t")
rds2loom.columns = ['rdscellname','loomcellname']
rds2loom = rds2loom[rds2loom.loomcellname.isin(CellID)]
umap = pd.read_csv(coord2d,header= 0,sep="\t")
umap.columns = ['rdscellname','UMAP_1','UMAP_2']
rds2loom = pd.merge(rds2loom,umap,on = "rdscellname")

cell_id_to_index = {cell_id: idx for idx, cell_id in enumerate(CellID)}
index = [cell_id_to_index.get(m, None) for m in rds2loom.loomcellname]
index = sorted(index) # loom要求index必须是递增的
index = np.array(index)
index_gene = np.array({idx for idx, gene_id in enumerate(Genename)})

with loompy.connect(singleloom) as ds:
  CellID = {'CellID':ds.ca.CellID[index].astype(str).tolist()}
  Gene = {'Gene':ds.ra.Gene[index_gene].astype(str).tolist()}
  ambiguous_full = ds.layers['ambiguous'][:, index][index_gene, :]
  spliced_full = ds.layers['spliced'][:, index][index_gene, :]
  unspliced_full  = ds.layers['unspliced'][:, index][index_gene, :]

# 创建obs层
obs_df = pd.DataFrame(index=CellID["CellID"])
obs_df.index.name = "CellID"

# 创建var层
var_df = pd.DataFrame(index=Gene["Gene"])
var_df.index.name = "Gene"

adata = ad.AnnData(
  X=csr_matrix(spliced_full.T),  # 填充 spliced 层数据
  obs=obs_df,  # 细胞名
  var=var_df,  # 基因名
  layers={
    "ambiguous": csr_matrix(ambiguous_full.T),
    "unspliced": csr_matrix(unspliced_full.T)
  }
)

# 替换为单细胞细胞名
adata = adata[adata.obs_names.isin(rds2loom.loomcellname)]
df1 = rds2loom.set_index("loomcellname")
df1.index.name = None
df1 = df1.loc[adata.obs_names.tolist(),:]
adata.obs_names = df1.loc[:,"rdscellname"].tolist()
adata.obs["clusters"] = df1.clusters.astype(str).tolist()
adata.obsm['X_umap'] = np.array(df1.loc[:,["UMAP_1","UMAP_2"]])

# 添加clusters
celllist = pd.read_csv(Celllist,sep= "\t",header = 0)
celllist = celllist.rename(columns = {str(celllist.columns[0]):"CellID"})
celllist = celllist.drop_duplicates(subset='CellID', keep='first', inplace=False)
mask = adata.obs_names.isin(celllist.iloc[:,0])
adata = adata[mask]
celllist = celllist.iloc[:,0:2]
celllist.columns = ['CellID', 'clusters']
celllist_indexed = celllist.set_index('CellID')
adata.obs['clusters'] = celllist_indexed.loc[adata.obs.index, 'clusters'].astype(str)

################################# scvelo分析代码

adata.var_names_make_unique()
scv.pp.remove_duplicate_cells(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts1, n_top_genes=n_top_genes1)
scv.pp.moments(adata, n_pcs=n_pcs1, n_neighbors=n_neighbors1)
scv.tl.umap(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata) 
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='clusters')
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)

