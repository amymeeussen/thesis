#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=dem_den_tax
#SBATCH --error=/home/aparsons/qiime/log/core.log
#SBATCH --output=/home/aparsons/qiime/log/core.log

echo "Activating qiime environment"
source activate qiime2-2022.8

#echo "Starting mouth core" && \
#date && \
#qiime feature-table core-features --i-table ~/qiime/ANCOM/mouth-table.qza --o-visualization ~/qiime/core/mouth-core.qzv && \
#qiime feature-table core-features --i-table ~/qiime/ANCOM/foot-table.qza --o-visualization ~/qiime/core/foot-core.qzv && \
#qiime feature-table core-features --i-table ~/qiime/ANCOM/cloaca-table.qza --o-visualization ~/qiime/core/cloaca-core.qzv && \
#\
#qiime metadata tabulate --m-input-file ~/qiime/core/mouth-core-features-80.tsv --o-visualization mouth-core-features-80.qzv && \
#qiime feature-table filter-features --i-table ~/qiime/ANCOM/mouth-table.qza --o-filtered-table ~/qiime/core/mouth-table-core-80.qza --m-metadata-file ~/qiime/core/mouth-core-features-80.tsv && \
#echo "starting barplot" && \
#qiime taxa barplot --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --i-table ~/qiime/core/mouth-table-core-80.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core/barplot-mouth-core-80.qzv && \
#\
#qiime metadata tabulate --m-input-file ~/qiime/core/foot-core-features-80.tsv --o-visualization foot-core-features-80.qzv && \
#qiime feature-table filter-features --i-table ~/qiime/ANCOM/foot-table.qza --o-filtered-table ~/qiime/core/foot-table-core-80.qza --m-metadata-file ~/qiime/core/foot-core-features-80.tsv && \
#echo "starting barplot" && \
#qiime taxa barplot --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --i-table ~/qiime/core/foot-table-core-80.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core/barplot-foot-core-80.qzv && \
#\
qiime metadata tabulate --m-input-file ~/qiime/core/cloaca-core-features-80.tsv --o-visualization cloaca-core-features-80.qzv && \
qiime feature-table filter-features --i-table ~/qiime/ANCOM/cloaca-table.qza --o-filtered-table ~/qiime/core/cloaca-table-core-80.qza --m-metadata-file ~/qiime/core/cloaca-core-features-80.tsv && \
echo "starting barplot" && \
qiime taxa barplot --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --i-table ~/qiime/core/cloaca-table-core-80.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/core/barplot-cloaca-core-80.qzv && \
date
