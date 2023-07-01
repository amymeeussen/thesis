#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=dem_den_tax
#SBATCH --error=/home/aparsons/qiime/log/whole3.log
#SBATCH --output=/home/aparsons/qiime/log/whole3.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8


#echo "Starting demultiplex" && \
#rm -rf ~/qiime/demux/ && \
#mkdir -p ~/qiime/demux && \
#date && \
# negative controls and birds 20, 60, 66 and 46 are filtered out of the metadata
#qiime demux emp-paired --m-barcodes-file ~/thesis/metadata/metadata.tsv --m-barcodes-column barcodes --i-seqs ~/qiime/lab_data/emp-paired-end-sequences.qza --o-per-sample-sequences ~/qiime/demux/demux.qza --o-error-correction-details ~/qiime/demux/demux-details.qza --p-no-golay-error-correction && \
#qiime demux summarize --i-data ~/qiime/demux/demux.qza --o-visualization ~/qiime/demux/demux.qzv && \
#echo "finished demultiplexing" && \
#\
#echo "Starting denoising" && \
#rm -rf ~/qiime/denoising/ && \
#mkdir -p ~/qiime/denoising && \
#date && \
#qiime dada2 denoise-paired --i-demultiplexed-seqs ~/qiime/demux/demux.qza --p-n-threads 24 --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 160 --p-trunc-len-r 160 --o-table ~/qiime/denoising/dada2-table.qza --o-representative-sequences ~/qiime/denoising/rep-seqs-dada2.qza --o-denoising-stats ~/qiime/denoising/denoising-stats.qza && \
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
#echo "start filter seq-rep" && \
# seq=reps: filter out Mitochondria, Chloroplast, Eukaryota/ include all taxa with phylum labels/ decontam table not used
#qiime taxa filter-seqs --i-sequences ~/qiime/denoising/rep-seqs-dada2.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include p__  --p-exclude Mitochondria,Chloroplast --o-filtered-sequences ~/qiime/denoising/sequences-filtered-noM-noC.qza && \
#qiime taxa filter-seqs --i-sequences ~/qiime/denoising/sequences-filtered-noM-noC.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include p__  --p-exclude "d_Eukaryota" --o-filtered-sequences ~/qiime/denoising/sequences-filtered.qza && \
#echo "finish filter seq-rep" && \
#\
#echo "start filter taxa table" && \
# taxa table: filter out Mitochondria, Chloroplast, Eukaryota/ include all taxa with phylum labels/ decontam table not used
#qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include p__ --p-exclude Mitochondria,Chloroplast --o-filtered-table ~/qiime/denoising/taxa-table-filtered-noM-noC.qza && \
#qiime taxa filter-table --i-table ~/qiime/denoising/taxa-table-filtered-noM-noC.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-exclude "d_Eukaryota" --o-filtered-table ~/qiime/denoising/feature-table-filtered.qza && \
#echo "finish filter taxa table" && \
#\
#echo "start taxa barplot" && \
#qiime taxa barplot --i-table ~/qiime/denoising/feature-table-filtered.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/denoising/taxa-bar-plot.qzv && \
#echo "finished taxa barplot" && \
#\
#echo "Starting PhylogeneticTree" && \
#rm -rf ~/qiime/PhylogeneticTree/ && \
#mkdir -p ~/qiime/PhylogeneticTree && \
#time qiime fragment-insertion sepp --i-representative-sequences ~/qiime/denoising/sequences-filtered.qza --i-reference-database ~/qiime/train_classifier/sepp-refs-silva-128.qza --o-tree ~/qiime/PhylogeneticTree/SILVAtree.qza --o-placements ~/qiime/PhylogeneticTree/SILVAtree_placements.qza --p-threads 20 && \
#qiime fragment-insertion filter-features --i-table ~/qiime/denoising/feature-table-filtered.qza --i-tree ~/qiime/PhylogeneticTree/SILVAtree.qza --o-filtered-table ~/qiime/PhylogeneticTree/SILVA-used-fragments.qza --o-removed-table ~/qiime/PhylogeneticTree/SILVA-table-removed-fragments.qza --verbose && \
#\
#echo "start rarify" && \
#rm -rf ~/qiime/diversity/ && \
#mkdir -p ~/qiime/diversity/ && \
#qiime diversity alpha-rarefaction --i-table ~/qiime/denoising/feature-table-filtered.qza --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-max-depth 60000 --p-steps 30 --p-iterations 20 --o-visualization ~/qiime/diversity/alpha-rarefaction-rooted.qzv && \
#\
echo "Starting alpha diversity" && \
rm -rf ~/qiime/core-metric-results/ && \
qiime diversity core-metrics-phylogenetic --i-phylogeny ~/qiime/PhylogeneticTree/SILVAtree.qza --i-table ~/qiime/denoising/feature-table-filtered.qza --p-sampling-depth 15000 --m-metadata-file ~/thesis/metadata/metadata.tsv --output-dir ~/qiime/core-metric-results/ && \ 
qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/faith-pd-group-significance.qzv && \
qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results/evenness_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/evenness-group-significance.qzv && \
qiime diversity alpha-correlation --i-alpha-diversity ~/qiime/core-metric-results/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/alpha-correlation-test.qzv && \
\
echo "start beta diversity" && \
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-Area-group-significance.qzv && \
qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column type --p-pairwise --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-sampletype-group-significance.qzv && \
\
qiime emperor plot --i-pcoa ~/qiime/core-metric-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-emperor.qzv && \
qiime emperor plot --i-pcoa ~/qiime/core-metric-results/bray_curtis_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core-metric-results/bray-curtis-emperor.qzv && \
\
echo "start mouth differential abundance" && \
rm -rf ~/qiime/ANCOM/ && \
mkdir -p ~/qiime/ANCOM && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where '[type]='"'"'M'"'"'' --o-filtered-table ~/qiime/ANCOM/mouth-table.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/ancom-mouth-Area.qzv && \
\
echo "start foot differential abundance" && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where '[type]='"'"'F'"'"'' --o-filtered-table ~/qiime/ANCOM/foot-table.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/foot-table.qza --o-composition-table ~/qiime/ANCOM/comp-foot-table.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-foot-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/ancom-foot-Area.qzv && \
\
echo "start cloaca differential abundance" && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where '[type]='"'"'C'"'"'' --o-filtered-table ~/qiime/ANCOM/cloaca-table.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/ancom-cloaca-Area.qzv && \
\
echo "start mouth differential abundance at genus level" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/mouth-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/mouth-table-l6.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table-l6.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l6-ancom-mouth-Area.qzv && \
\
echo "start foot differential abundance at genus level" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/foot-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/foot-table-l6.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/foot-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-foot-table-l6.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-foot-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l6-ancom-foot-Area.qzv && \
\
echo "start cloaca differential abundance at genus level" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/cloaca-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/cloaca-table-l6.qza && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l6-ancom-cloaca-Area.qzv && \
\
echo "Start identifying core microbiome" && \
rm -rf ~/qiime/core/ && \
mkdir -p ~/qiime/core && \
date && \
qiime feature-table core-features --i-table ~/qiime/denoising/feature-table-filtered.qza --o-visualization ~/qiime/core/core.qzv && \
date
