#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=decontam
#SBATCH --error=/home/aparsons/qiime/log/decontam.log
#SBATCH --output=/home/aparsons/qiime/log/decontam.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8
#echo "starting clean" && \
#qiime taxa filter-table --i-table ~/qiime/denoising/decontam_table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-exclude 'p_; ,Mitochondria,Chloroplast,Eukaryota' --o-filtered-table ~/qiime/denoising/decontam_table_filtered.qza && \
#echo "finish clean" && \
#\
#echo "start sequence reclassification" && \
#qiime feature-table filter-seqs --i-data ~/qiime/denoising/rep-seqs-dada2.qza --i-table ~/qiime/denoising/decontam_table_filtered.qza --o-filtered-data ~/qiime/denoising/decontam_sequences_filtered.qza && \
#echo "finish sequence reclassification" && \
#\
#echo "start taxa barplot" && \
#qiime taxa barplot --i-table ~/qiime/denoising/decontam_table_filtered.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/denoising/decontam-taxa-bar-plot-filtered.qzv && \
#echo "finished taxa barplot" && \
#\
#echo "Starting removing unwanted" && \
#time qiime fragment-insertion sepp --i-representative-sequences ~/qiime/denoising/decontam_sequences_filtered.qza --i-reference-database ~/qiime/train_classifier/sepp-refs-silva-128.qza --o-tree ~/qiime/PhylogeneticTree/decontam_SILVAtree.qza --o-placements ~/qiime/PhylogeneticTree/decontam_SILVAtree_placements.qza$
#qiime fragment-insertion filter-features --i-table ~/qiime/denoising/decontam_table_filtered.qza --i-tree ~/qiime/PhylogeneticTree/decontam_SILVAtree.qza --o-filtered-table ~/qiime/PhylogeneticTree/decontam_SILVA-filtered-table.qza --o-removed-table ~/qiime/PhylogeneticTree/decontam_SILVA_removed_table.qza$
#/
#echo "start rarifying"
#qiime diversity alpha-rarefaction --i-table ~/qiime/denoising/decontam_table_filtered.qza --i-phylogeny ~/qiime/PhylogeneticTree/decontam_SILVAtree.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-max-depth 29652 --p-steps 30 --p-iterations 20 --o-visualization ~/qiime/diversity/decontam_alpha-rarefaction-rooted.qzv && \
#/
#echo "Starting alpha diversity" && \
#rm -rf ~/qiime/decontam-core-metric-results/ && \
#qiime diversity core-metrics-phylogenetic --i-phylogeny ~/qiime/PhylogeneticTree/decontam_SILVAtree.qza --i-table ~/qiime/PhylogeneticTree/decontam_SILVA-filtered-table.qza --p-sampling-depth 29652 --m-metadata-file ~/thesis/metadata/metadata.tsv --output-dir ~/qiime/decontam-core-metric-results/ && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/decontam-core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/decontam-core-metric-results/decontam-faith-pd-group-significance.qzv && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/decontam-core-metric-results/evenness_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/decontam-core-metric-results/decontam-evenness-group-significance.qzv && \
#qiime diversity alpha-correlation --i-alpha-diversity ~/qiime/decontam-core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/decontam-core-metric-results/decontam-alpha-correlation-test.qzv && \
#\
echo "start beta analysis" && \
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/decontam-core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/decontam-core-metric-results//decontam-unweighted-unifrac-Area-group-significance.qzv && \
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column type --p-pairwise --o-visualization ~/qiime/decontam-core-metric-results/decontam-unweighted-unifrac-sampletype-group-significance.qzv && \
\
qiime emperor plot --i-pcoa ~/qiime/decontam-core-metric-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/decontam-core-metric-results/decontam-unweighted-unifrac-emperor.qzv && \
qiime emperor plot --i-pcoa ~/qiime/decontam-core-metric-results/bray_curtis_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/decontam-core-metric-results/decontam-bray-curtis-emperor.qzv && \
date\


