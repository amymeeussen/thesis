#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=beta
#SBATCH --error=beta.log
#SBATCH --output=beta.log

date && \
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "Starting group significance" && \
qiime diversity beta-group-significance --i-distance-matrix core-metric-results/unweighted-unifrac-distance-matrix.qza --m-metadata-file metadata.tsv --m-metadata-column "Area" --o-visualization core-metric-results/unweighted-unifrac-body-site-significance.qzv --p-pairwise && \
echo "Starting beta group significance" && \
qiime diversity beta-group-significance --i-distance-matrix core-metric-results/unweighted-unifrac-distance-matrix.qza --m-metadata-file metadata.tsv --m-metadata-column "Area" --o-visualization core-metric-results/unweighted-unifrac-subject-group-significance.qzv --p-pairwise && \
date \







