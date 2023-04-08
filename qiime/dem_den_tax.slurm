#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=dem_den_tax
#SBATCH --error=/home/aparsons/qiime/log/dem_den_tax.log
#SBATCH --output=/home/aparsons/qiime/log/dem_den_tax.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8


#echo "Starting demultiplexing" && \
#rm -rf ~/qiime/demux/ && \
#mkdir -p ~/qiime/demux && \
#date && \
#qiime demux emp-paired --m-barcodes-file ~/thesis/metadata/metadata.tsv --m-barcodes-column barcodes --i-seqs ~/qiime/lab_data/emp-paired-end-sequences.qza --o-per-sample-sequences ~/qiime/demux/demux.qza --o-error-correction-details ~/qiime/demux/demux-details.qza --p-no-golay-error-correction && \
#qiime demux summarize --i-data ~/qiime/demux/demux.qza --o-visualization ~/qiime/demux/demux.qzv && \
#echo "finished demultiplexing" && \
#\
#echo "Starting denoising" && \
#rm -rf ~/qiime/denoising/ && \
#mkdir -p ~/qiime/denoising && \
#date && \
#qiime dada2 denoise-paired --i-demultiplexed-seqs ~/qiime/demux/demux.qza --p-n-threads 24 --p-trim-left-f 13 --p-trim-left-r 7 --p-trunc-len-f 170 --p-trunc-len-r 169 --o-table ~/qiime/denoising/dada2-table.qza --o-representative-sequences ~/qiime/denoising/rep-seqs-dada2.qza --o-denoising-stats ~/qiime/denoising/denoising-stats.qza && \
#echo "Finished denoising" && \
#echo "Starting denoisingstats" && \
#date && \
#qiime feature-table summarize --i-table ~/qiime/denoising/dada2-table.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/denoising/dada2-table.qzv && \
#qiime feature-table tabulate-seqs --i-data ~/qiime/denoising/rep-seqs-dada2.qza --o-visualization ~/qiime/denoising/rep-seqs.qzv && \
#qiime metadata tabulate --m-input-file ~/qiime/denoising/denoising-stats.qza --o-visualization ~/qiime/denoising/denoising-stats.qzv && \
#echo "Finished denoisingstats" && \
#\
#echo "starting taxonomy assignment" && \
#rm -rf ~/qiime/taxonomy/ && \
#mkdir -p ~/qiime/taxonomy && \
#qiime feature-classifier classify-sklearn --i-classifier ~/qiime/train_classifier/classifier.qza --i-reads ~/qiime/denoising/rep-seqs-dada2.qza --verbose --o-classification ~/qiime/taxonomy/taxonomy.qza && \
#echo "finish taxonomy assignment" && \
#date && \
#echo "start visualization" && \
#qiime metadata tabulate --m-input-file ~/qiime/taxonomy/taxonomy.qza --o-visualization ~/qiime/taxonomy/taxonomy.qzv && \
#echo "finished visualization" && \
#date && \
#echo "start barplot" && \
#qiime taxa barplot --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/taxonomy/taxa-bar-plots.qzv && \
#echo "finished barplot" && \
#\
#echo "starting clean" && \
#qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza/ --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-exclude 'p_; ,Mitochondria,Chloroplast,Eukaryota' --o-filtered-table ~/qiime/denoising/table-filtered.qza && \
#echo "finish clean" \
#&& \
#echo "start reclassification" && \
#qiime feature-table filter-seqs --i-data ~/qiime/denoising/rep-seqs-dada2.qza --i-table ~/qiime/denoising/table-filtered.qza --o-filtered-data ~/qiime/denoising/sequences-filtered.qza && \
#echo "finish reclassification" \
#&& \
#echo "start taxa barplot" && \
#qiime taxa barplot --i-table ~/qiime/denoising/table-filtered.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/denoising/taxa-bar-plot-filtered.qzv && \
#echo "finished taxa barplot" && \
#\
#echo "Starting removing unwanted" && \
#time qiime fragment-insertion sepp --i-representative-sequences ~/qiime/denoising/sequences-filtered.qza --i-reference-database ~/qiime/train_classifier/sepp-refs-silva-128.qza --o-tree ~/qiime/PhylogeneticTree/SILVAtree.qza --o-placements ~/qiime/PhylogeneticTree/SILVAtree_placements.qza --p-threads 20
#qiime fragment-insertion filter-features --i-table ~/qiime/denoising/table-filtered.qza --i-tree ~/qiime/PhylogeneticTree/SILVAtree.qza --o-filtered-table ~/qiime/PhylogeneticTree/SILVA-filtered-table.qza --o-removed-table ~/qiime/PhylogeneticTree/removed_table.qza --verbose && \
#\
#echo "start rarify"
#qiime diversity alpha-rarefaction --i-table ~/qiime/denoising/table-filtered.qza --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-max-depth 29652 --p-steps 30 --p-iterations 20 --o-visualization ~/qiime/diversity/alpha-rarefaction-rooted.qzv
#\
#echo "Starting alpha diversity" && \
#rm -rf ~/qiime/core-metric-results/ && \
#qiime diversity core-metrics-phylogenetic --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --i-table ~/qiime/PhylogeneticTree/SILVA-filtered-table.qza --p-sampling-depth 29652 --m-metadata-file ~/thesis/metadata/metadata.tsv --output-dir ~/qiime/core-metric-results/ 
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/decontam-faith-pd-group-significance.qzv && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/evenness_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/decontam-evenness-group-significance.qzv && \
#qiime diversity alpha-correlation --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/decontam-alpha-correlation-test.qzv && \
#\
echo "start beta diversity"
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-Area-group-significance.qzv && \
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column type --p-pairwise --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-sampletype-group-significance.qzv && \
\
qiime emperor plot --i-pcoa ~/qiime/core-metric-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-emperor.qzv && \
qiime emperor plot --i-pcoa ~/qiime/core-metric-results/bray_curtis_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/bray-curtis-emperor.qzv && \
date\
