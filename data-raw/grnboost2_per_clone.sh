#!/bin/bash

mkdir -p $3

for filepath in $1/*
do
  for i in 1 2 3 4 5 6 7 8 9 10
  do
    sbatch --job-name="GRNBoost2_"$(basename "$filepath" .csv) --exclude gorilla4 --wrap=". /opt/miniconda3/etc/profile.d/conda.sh
                   conda activate GRNBoost2
                   python data-raw/arboreto_with_multiprocessing.py \
                          $filepath \
                          $2 \
                          --method grnboost2 \
                          --output $3"/"$(basename "$filepath" .csv)"_"$i".tsv" \
                          --num_workers 8 \
                          --seed $i"
  done
done

# to run the script from commandline:
# sbatch grnboost2_per_clone.sh "/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.network_inference/input/count_matrices" "/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.network_inference/input/regulators.txt" "/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.network_inference/output"
