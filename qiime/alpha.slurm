#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=alpha
#SBATCH --error=alpha.log
#SBATCH --output=alpha.log

date \
&& \
echo "Activating qiime environemnt" \
&& \
source activate qiime2-2022.8 \
&& \
##echo "Starting alpha diversity" \
##&& \
##qiime diversity core-metrics-phylogenetic --i-phylogeny PhylogeneticTree/SILVAtree.qza --i-table PhylogeneticTree/SILVA-filtered-table.qza --p-sampling-depth 20000 --m-metadata-file metadata.tsv --o-rarefied-table core-metric-results/rarefied_table.qza --o-observed-features-vector core-metric-results/observed-otus-vector.qza --o-shannon-vector core-metric-results/shannon-vector.qza --o-jaccard-distance-matrix core-metric-results/jaccard-distance-matrix.qza --o-bray-curtis-distance-matrix core-metric-results/bray-curtis-distance-matrix.qza --o-jaccard-pcoa-results core-metric-results/jaccard-pcoa-results.qza --o-bray-curtis-pcoa-results core-metric-results/bray-curtis-pcoa-results.qza --o-jaccard-emperor core-metric-results/jaccard-emperor.qzv --o-bray-curtis-emperor core-metric-results/bray-curtis-emperor.qzv --o-faith-pd-vector core-metric-results/faith-pd-vector.qza --o-evenness-vector core-metric-results/evenness-vector.qza --o-unweighted-unifrac-distance-matrix core-metric-results/unweighted-unifrac-distance-matrix.qza --o-weighted-unifrac-distance-matrix core-metric-results/weighted-unifrac-distance-matrix.qza --o-unweighted-unifrac-pcoa-results core-metric-results/unweighted-unifrac-pcoa-results.qza --o-unweighted-unifrac-emperor core-metric-results/unweighted-unifrac-emperor.qzv --o-weighted-unifrac-emperor core-metric-results/weighted-unifrac-empeor.qzv --o-weighted-unifrac-pcoa-results core-metric-results/weighted-unifrac-pcoa-results.qza \
##&& \
##echo "Finishing alpha diversity" \
##&& \
##echo "Starting faith_pd_vector" \
##&& \
##qiime diversity alpha-group-significance --i-alpha-diversity core-metric-results/faith-pd-vector.qza --m-metadata-file metadata.tsv --o-visualization core-metric-results/faith-pd-group-significance.qzv \
##&& \
##qiime diversity alpha-group-significance --i-alpha-diversity core-metric-results/evenness-vector.qza --m-metadata-file metadata.tsv --o-visualization core-metric-results/evenness-group-significance.qzv \
##&& \
qiime diversity alpha-correlation --i-alpha-diversity core-metric-results/faith-pd-vector.qza --m-metadata-file metadata.tsv --o-visualization alpha-correlation-test.qza \
&& \
date \


