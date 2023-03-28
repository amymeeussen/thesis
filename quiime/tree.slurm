#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=tree
#SBATCH --error=tree.log
#SBATCH --output=tree.log

date \
&& \
echo "Activating qiime environment" \
&& \
source activate qiime2-2022.8 \
&& \
echo "Starting removing unwanted" && \
time qiime fragment-insertion sepp --i-representative-sequences filtered-sequences.qza --i-reference-database sepp-refs-silva-128.qza --o-tree PhylogeneticTree/SILVAtree.qza --o-placements SILVAtree_placements.qza --p-threads 10 && \
qiime fragment-insertion filter-features --i-table table-filtered.qza --i-tree PhylogeneticTree/SILVAtree.qza --o-filtered-table PhylogeneticTree/SILVA-filtered-table.qza --o-removed-table PhylogeneticTree/removed_table.qza --verbose && \
date \

