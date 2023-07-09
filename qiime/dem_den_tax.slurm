#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=dem_den_tax_ss
#SBATCH --error=/home/aparsons/qiime/log/whole3_ss.log
#SBATCH --output=/home/aparsons/qiime/log/whole3_ss.log

export SUFFIX="_ss"
echo "Using suffix " ${SUFFIX}

echo "Activating qiime environemnt"
source activate qiime2-2022.8

#echo "Starting demultiplex" && \
#rm -rf ~/qiime/demux${SUFFIX}/ && \
#mkdir -p ~/qiime/demux${SUFFIX} && \
#date && \
# After looking at rarefication curves, I took birds 46 and 66 out of the dataset, as well as 15 other  birds to make equal sample sizes for ML and SF
#qiime demux emp-paired --m-barcodes-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-barcodes-column barcodes --i-seqs ~/qiime/lab_data/emp-paired-end-sequences.qza --o-per-sample-sequences ~/qiime/demux${SUFFIX}/demux.qza --o-error-correction-details ~/qiime/demux${SUFFIX}/demux-details.qza --p-no-golay-error-correction && \
#qiime demux summarize --i-data ~/qiime/demux${SUFFIX}/demux.qza --o-visualization ~/qiime/demux${SUFFIX}/demux${SUFFIX}.qzv && \
#echo "finished demultiplexing" && \
#\
#echo "Starting denoising" && \
#rm -rf ~/qiime/denoising${SUFFIX}/ && \
#mkdir -p ~/qiime/denoising${SUFFIX} && \
#date && \
#qiime dada2 denoise-paired --i-demultiplexed-seqs ~/qiime/demux${SUFFIX}/demux.qza --p-n-threads 24 --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 160 --p-trunc-len-r 160 --o-table ~/qiime/denoising${SUFFIX}/dada2-table.qza --o-representative-sequences ~/qiime/denoising${SUFFIX}/rep-seqs-dada2.qza --o-denoising-stats ~/qiime/denoising${SUFFIX}/denoising-stats.qza && \
#echo "Finished denoising" && \
#echo "Starting denoisingstats" && \
#date && \
#qiime feature-table summarize --i-table ~/qiime/denoising${SUFFIX}/dada2-table.qza --m-sample-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/denoising${SUFFIX}/dada2-table${SUFFIX}.qzv && \
#qiime feature-table tabulate-seqs --i-data ~/qiime/denoising${SUFFIX}/rep-seqs-dada2.qza --o-visualization ~/qiime/denoising${SUFFIX}/rep-seqs${SUFFIX}.qzv && \
#qiime metadata tabulate --m-input-file ~/qiime/denoising${SUFFIX}/denoising-stats.qza --o-visualization ~/qiime/denoising${SUFFIX}/denoising-stats${SUFFIX}.qzv && \
#echo "Finished denoisingstats" && \
#\
#echo "starting taxonomy assignment" && \
#rm -rf ~/qiime/taxonomy${SUFFIX}/ && \
#mkdir -p ~/qiime/taxonomy${SUFFIX} && \
#qiime feature-classifier classify-sklearn --i-classifier ~/qiime/train_classifier/classifier.qza --i-reads ~/qiime/denoising${SUFFIX}/rep-seqs-dada2.qza --verbose --o-classification ~/qiime/taxonomy${SUFFIX}/taxonomy.qza && \
#echo "finish taxonomy assignment" && \
#date && \
#echo "start visualization" && \
#qiime metadata tabulate --m-input-file ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --o-visualization ~/qiime/taxonomy${SUFFIX}/taxonomy${SUFFIX}.qzv && \
#echo "finished visualization" && \
#date && \
#echo "start filter seq-rep" && \
# seq=reps: filter out Mitochondria, Chloroplast, Eukaryota/ include all taxa with phylum labels/ decontam table not used
#qiime taxa filter-seqs --i-sequences ~/qiime/denoising${SUFFIX}/rep-seqs-dada2.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-include p__  --p-exclude Mitochondria,Chloroplast --o-filtered-sequences ~/qiime/denoising${SUFFIX}/sequences-filtered-noM-noC.qza && \
#qiime taxa filter-seqs --i-sequences ~/qiime/denoising${SUFFIX}/sequences-filtered-noM-noC.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-include p__  --p-exclude "d_Eukaryota" --o-filtered-sequences ~/qiime/denoising${SUFFIX}/sequences-filtered.qza && \
#echo "finish filter seq-rep" && \
#\
#echo "start filter taxa table" && \
#taxa table: filter out Mitochondria, Chloroplast, Eukaryota/ include all taxa with phylum labels/ decontam table not used
#qiime taxa filter-table --i-table ~/qiime/denoising${SUFFIX}/dada2-table.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-include p__ --p-exclude Mitochondria,Chloroplast --o-filtered-table ~/qiime/denoising${SUFFIX}/taxa-table-filtered-noM-noC.qza && \
#qiime taxa filter-table --i-table ~/qiime/denoising${SUFFIX}/taxa-table-filtered-noM-noC.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-exclude "d_Eukaryota" --o-filtered-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza && \
#echo "finish filter taxa table" && \
#\
#echo "start taxa barplot" && \
#qiime taxa barplot --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/denoising${SUFFIX}/taxa-bar-plot${SUFFIX}.qzv && \
#echo "finished taxa barplot" && \
#\
#echo "Starting PhylogeneticTree" && \
#rm -rf ~/qiime/PhylogeneticTree${SUFFIX}/ && \
#mkdir -p ~/qiime/PhylogeneticTree${SUFFIX} && \
#time qiime fragment-insertion sepp --i-representative-sequences ~/qiime/denoising${SUFFIX}/sequences-filtered.qza --i-reference-database ~/qiime/train_classifier/sepp-refs-silva-128.qza --o-tree ~/qiime/PhylogeneticTree${SUFFIX}/SILVAtree.qza --o-placements ~/qiime/PhylogeneticTree${SUFFIX}/SILVAtree_placements.qza --p-threads 20 && \
#qiime fragment-insertion filter-features --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --i-tree ~/qiime/PhylogeneticTree${SUFFIX}/SILVAtree.qza --o-filtered-table ~/qiime/PhylogeneticTree${SUFFIX}/SILVA-used-fragments.qza --o-removed-table ~/qiime/PhylogeneticTree${SUFFIX}/SILVA-table-removed-fragments.qza --verbose && \
#\
#echo "start rarify" && \
#rm -rf ~/qiime/diversity${SUFFIX}/ && \
#mkdir -p ~/qiime/diversity${SUFFIX}/ && \
#qiime diversity alpha-rarefaction --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --i-phylogeny ~/qiime/PhylogeneticTree${SUFFIX}/SILVAtree.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-max-depth 60000 --p-steps 30 --p-iterations 20 --o-visualization ~/qiime/diversity${SUFFIX}/alpha-rarefaction-rooted${SUFFIX}.qzv && \
#\
#echo "Starting alpha diversity" && \
#rm -rf ~/qiime/core-metric-results${SUFFIX}/ && \
#qiime diversity core-metrics-phylogenetic --i-phylogeny ~/qiime/PhylogeneticTree${SUFFIX}/SILVAtree.qza --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --p-sampling-depth 10000 --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --output-dir ~/qiime/core-metric-results${SUFFIX}/ && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results${SUFFIX}/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/core-metric-results${SUFFIX}/faith-pd-group-significance${SUFFIX}.qzv && \
#qiime diversity alpha-group-significance --i-alpha-diversity ~/qiime/core-metric-results${SUFFIX}/evenness_vector.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/core-metric-results${SUFFIX}/evenness-group-significance${SUFFIX}.qzv && \
#qiime diversity alpha-correlation --i-alpha-diversity ~/qiime/core-metric-results${SUFFIX}/faith_pd_vector.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/core-metric-results${SUFFIX}/alpha-correlation-test${SUFFIX}.qzv && \
#\
#echo "start beta diversity" && \
#qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results${SUFFIX}/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/core-metric-results${SUFFIX}/unweighted-unifrac-Area-group-significance${SUFFIX}.qzv && \
#qiime diversity beta-group-significance --i-distance-matrix ~/qiime/core-metric-results${SUFFIX}/unweighted_unifrac_distance_matrix.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column type --p-pairwise --o-visualization ~/qiime/core-metric-results${SUFFIX}/unweighted-unifrac-sampletype-group-significance${SUFFIX}.qzv && \
#\
#qiime emperor plot --i-pcoa ~/qiime/core-metric-results${SUFFIX}/unweighted_unifrac_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/core-metric-results/unweighted-unifrac-emperor${SUFFIX}.qzv && \
#qiime emperor plot --i-pcoa ~/qiime/core-metric-results${SUFFIX}/bray_curtis_pcoa_results.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --o-visualization ~/qiime/core-metric-results/bray-curtis-emperor${SUFFIX}.qzv && \
#\
#echo "start mouth differential abundance" && \
#rm -rf ~/qiime/ANCOM${SUFFIX}/ && \
#mkdir -p ~/qiime/ANCOM${SUFFIX} && \
#qiime feature-table filter-samples --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[type]='"'"'M'"'"'' --o-filtered-table ~/qiime/ANCOM${SUFFIX}/mouth-table.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/mouth-table.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-mouth-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-mouth-table.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/ancom-mouth-Area${SUFFIX}.qzv && \
#\
#echo "start foot differential abundance" && \
#qiime feature-table filter-samples --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[type]='"'"'F'"'"'' --o-filtered-table ~/qiime/ANCOM${SUFFIX}/foot-table.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/foot-table.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-foot-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-foot-table.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/ancom-foot-Area${SUFFIX}.qzv && \
#\
#echo "start cloaca differential abundance" && \
#qiime feature-table filter-samples --i-table ~/qiime/denoising${SUFFIX}/feature-table-filtered.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[type]='"'"'C'"'"'' --o-filtered-table ~/qiime/ANCOM${SUFFIX}/cloaca-table.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-cloaca-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-cloaca-table.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/ancom-cloaca-Area${SUFFIX}.qzv && \
#\
#echo "start mouth differential abundance at genus level" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM${SUFFIX}/mouth-table.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM${SUFFIX}/mouth-table-l6.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/mouth-table-l6.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-mouth-table-l6.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/l6-ancom-mouth-Area${SUFFIX}.qzv && \
#\
#echo "start foot differential abundance at genus level" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM${SUFFIX}/foot-table.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-foot-table-l6.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-foot-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/l6-ancom-foot-Area${SUFFIX}.qzv && \
#\
#echo "start cloaca differential abundance at genus level" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table.qza --i-taxonomy ~/qiime/taxonomy${SUFFIX}/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6.qza && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6.qza --o-composition-table ~/qiime/ANCOM${SUFFIX}/comp-cloaca-table-l6.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM${SUFFIX}/comp-cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM${SUFFIX}/l6-ancom-cloaca-Area${SUFFIX}.qzv && \
#\
#echo "Start identifying core microbiome" && \
#rm -rf ~/qiime/core${SUFFIX}/ && \
#mkdir -p ~/qiime/core${SUFFIX} && \
#\
# Mouth core from ML
#qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'ML'"'"'' --o-filtered-table ~/qiime/ANCOM_ss/mouth-table-l6_ML.qza && \
#qiime feature-table core-features --i-table ~/qiime/ANCOM_ss/mouth-table-l6_ML.qza  --o-visualization ~/qiime/core${SUFFIX}/mouth_ML_core${SUFFIX}.qzv && \
#\
# Mouth core from SF
#qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'SF'"'"'' --o-filtered-table ~/qiime/ANCOM_ss/mouth-table-l6_SF.qza && \
#qiime feature-table core-features --i-table ~/qiime/ANCOM_ss/mouth-table-l6_SF.qza  --o-visualization ~/qiime/core${SUFFIX}/mouth_SF_core${SUFFIX}.qzv && \
#\
# Foot core from ML
qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'ML'"'"'' --o-filtered-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6_ML.qza && \
qiime feature-table core-features --i-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6_ML.qza  --o-visualization ~/qiime/core${SUFFIX}/foot_ML_core${SUFFIX}.qzv && \
\
# Foot core from SF
qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'SF'"'"'' --o-filtered-table ~/qiime/ANCOM_ss/foot-table-l6_SF.qza && \
qiime feature-table core-features --i-table ~/qiime/ANCOM${SUFFIX}/foot-table-l6_SF.qza  --o-visualization ~/qiime/core${SUFFIX}/foot_SF_core${SUFFIX}.qzv && \
\
# Cloaca core from ML
qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'ML'"'"'' --o-filtered-table ~/qiime/ANCOM_ss/cloaca-table-l6_ML.qza && \
qiime feature-table core-features --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6_ML.qza  --o-visualization ~/qiime/core${SUFFIX}/cloaca_ML_core${SUFFIX}.qzv && \
\
# Cloaca core from SF
qiime feature-table filter-samples --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata${SUFFIX}.tsv --p-where '[Area]='"'"'SF'"'"'' --o-filtered-table ~/qiime/ANCOM_ss/cloaca-table-l6_SF.qza && \
qiime feature-table core-features --i-table ~/qiime/ANCOM${SUFFIX}/cloaca-table-l6_SF.qza  --o-visualization ~/qiime/core${SUFFIX}/cloaca_SF_core${SUFFIX}.qzv && \
date \
