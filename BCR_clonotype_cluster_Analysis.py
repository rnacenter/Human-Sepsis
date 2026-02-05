#scirpy版本为0.11.1

contigfile = "vdj_b.filtered_contig_annotations.csv" # 过滤后的contig.csv文件
receptor_type_cut = "multichain,ambiguous,TCR,no_IR" # 过滤条目组1
receptor_subtype_cut = "TRG+TRD,TRA+TRB,multichain,ambiguous,no_IR" # 过滤条目组2
chain_pairing_cut = "multichain,ambiguous,no_IR" # 过滤条目组3

receptor_arms="all" # 计算时用的V(D)J基因组成
dual_ir="primary_only" # 是否考虑双受体
within_group = "receptor_type"# 对clonotype_id鉴定的类型约束
RNAfile = "scRNA.h5ad" # 对应单细胞数据
cutoff = 10 # alignment算法的cutoff阈值

fdr_correction = True # 是否记录clonetype模块化分数的fdr值(否则记录p值)
random_state = 0 #重复计算模块化分数的seed
n_neighbors=30 #计算克隆型模块化分数使用的nneighbors
n_pcs = 50 #计算克隆型模块化分数使用的PCs

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt, cm as mpl_cm
import cycler
import os
import pickle

adata = ir.io.read_10x_vdj(contigfile)
#################1.质控
ir.tl.chain_qc(adata)
receptor_type_cut = receptor_type_cut.split(",")
receptor_subtype_cut = receptor_subtype_cut.split(",")
chain_pairing_cut = chain_pairing_cut.split(",")
adata = adata[~adata.obs["receptor_type"].isin(receptor_type_cut), :]
adata = adata[~adata.obs["receptor_subtype"].isin(receptor_subtype_cut), :]
adata = adata[~adata.obs["chain_pairing"].isin(chain_pairing_cut), :]

#################2.clonotype_cluster鉴定
ir.pp.ir_dist(adata,metric = 'alignment',cutoff = cutoff,sequence = 'aa')
ir.tl.define_clonotypes(
  adata,   
  sequence='aa',
  metric='alignment',
  receptor_arms=receptor_arms, 
  dual_ir=dual_ir,
  inplace = True,
  same_v_gene = False,
  within_group = within_group,
  distance_key = "ir_dist_aa_alignment"
)
ir.tl.clonotype_network(
  adata,
  sequence='aa', 
  metric='alignment',
  clonotype_key = "cc_aa_alignment"
)
#################3.合并单细胞数据
adata_scRNAseq = sc.read(RNAfile)
ir.pp.merge_with_ir(adata_scRNAseq,adata)
adata_scRNAseq.uns = adata.uns
adata = adata_scRNAseq
# 细胞数变化,重算布局
ir.tl.clonotype_network(
  adata, 
  sequence='aa', 
  metric='alignment',
  clonotype_key = "cc_aa_alignment"
)
ir.tl.clonal_expansion(adata,target_col = "cc_aa_alignment")

#############4.组间TCR信息比较
ir.tl.repertoire_overlap(adata,target_col = "cc_aa_alignment", groupby = "clusters", inplace=True,added_key = "repertoire_overlap") 
ir.tl.repertoire_overlap(adata,target_col = "cc_aa_alignment", groupby = "sample"
             , inplace=True,added_key = "repertoire_overlap_sample")

##############5.克隆型模块化分数 
sc.pp.neighbors(adata, n_neighbors=n_neighbors,n_pcs = n_pcs)
ir.tl.clonotype_modularity(adata,target_col="cc_aa_alignment",fdr_correction = fdr_correction,key_added = "clonotype_modularity",random_state = random_state) 

#############6.clonotype_size差异 
ir.tl.repertoire_overlap(adata, "clusters",target_col="cc_aa_alignment")
ir.tl.repertoire_overlap(adata, "sample",target_col="cc_aa_alignment")

####################7.保存
obs = adata.obs
obs.to_csv("BCRGeneAnno_cc_result.csv",sep = "\t")
pickle.dump(adata,open("BCRGeneAnno_cC_adata.dat","wb"))
