#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=mouth_alpha_beta
#SBATCH --error=/home/aparsons/qiime/log/mouth_alpha_beta.log
#SBATCH --output=/home/aparsons/qiime/log/mouth_alpha_beta.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8
echo "filter out mouth samples"
qiime feature-table filter-samples --i-table ~/qiime/denoising/table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where '[type]='"'"'M'"'"'' --o-filtered-table ~/qiime/denoising/mouth-table.qza && \
\
echo "start rarify"
qiime diversity alpha-rarefaction --i-table ~/qiime/denoising/mouth-table.qza --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-max-depth 60000 --p-steps 30 --p-iterations 20 --o-visualization ~/qiime/diversity/mouth_alpha-rarefaction-rooted.qzv
#\
#echo "Starting alpha diversity" && \
#rm -rf ~/qiime/core-metric-results/ &&
#qiime diversity core-metrics-phylogenetic --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --i-table ~/qiime/PhylogeneticTree/SILVA-filtered-table.qza --p-sampling-depth 29652 --m-metadata-file ~/thesis/metadata/metadata.tsv --output-dir$
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/faith-pd-group-significance.qzv && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/evenness_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/evenness-group-significance.qzv && \
#qiime diversity alpha-correlation --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/alpha-correlation-test.qzv && \
#\
#echo "start beta diversity"
#qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/core-metric-re$
#qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column type --p-pairwise --o-visualization ~/qiime/c$
#\
#qiime emperor plot --i-pcoa ~/qiime/core-metric-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-emperor.qzv && \
#qiime emperor plot --i-pcoa ~/qiime/core-metric-results/bray_curtis_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/bray-curtis-emperor.qzv && \
