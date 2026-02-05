#!/bin/bash
matrix="All.NormalizedCounts.csv"
tf="TFs.txt"
outfile="All.grn_output.tsv"
thread=1
#grnboost2或genie3
method="grnboost2"

pyscenic grn -o ${outfile} -t -m ${method} --num_workers ${thread} ${matrix} ${tf} --seed 2022
echo "grn finish"
pyscenic ctx All.grn_output.tsv "10kb_up_and_down_tss.feather" "500bp_up_and_100bp_down_tss.feather" -t -o "All.ctx_output.tsv" --num_workers 8 --nes_threshold 3 --auc_threshold 0.05 --rank_threshold 5000 --annotations_fname "Anno.txt" --expression_mtx_fname "All.NormalizedCounts.csv"
echo "ctx finish"
pyscenic aucell "All.NormalizedCounts.csv" "All.ctx_output.tsv" -t -o "All.aucell_output.tsv" --seed 2022 --num_workers 8
echo "aucell finish"