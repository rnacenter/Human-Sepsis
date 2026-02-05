#cellphone
cpdb_file_path = "cpdbpath/"                #cellphone数据库zip路径
meta_file_path = "celllist.txt"             #单细胞分析分群信息
counts_file_path = "args.countsfile"        #单细胞分析count表达矩阵(.txt)
out_path = "output_dir"                     #输出路径
precision = 3                               #计算结果精度
threshold = 0.1                             #表达细胞通讯基因的细胞占比阈值
iterations = 1000                           #计算p值迭代次数
thread = 4                                  #分析使用线程

import pandas as pd
import sys
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
os.chdir(out_path)
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 
    meta_file_path = meta_file_path,                 
    counts_file_path = counts_file_path,             
    counts_data = 'ensembl',                     
    microenvs_file_path = 'None',       
    iterations = iterations,                               
    threshold = threshold,                                 
    threads = thread,                                     
    debug_seed = 42,                                 
    result_precision = precision,                   
    pvalue = 0.05,                                  
    subsampling = False,                             
    subsampling_log = False,                         
    subsampling_num_pc = 100,                        
    subsampling_num_cells = 1000,                    
    separator = '_',                                 
    debug = False,                                   
    output_path = out_path,                          
    output_suffix = 'TestProject'                      
    )
