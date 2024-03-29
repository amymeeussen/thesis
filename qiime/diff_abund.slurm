#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=diff_abund_fam
#SBATCH --error=/home/aparsons/qiime/log/diff_abund_5.log
#SBATCH --output=/home/aparsons/qiime/log/diff_abund_5.log

# One of the assumptions of the differential abundance tests is that there is less that 25% difference btw two groups. Therefore, this file runs a differential abundance test on each of the sample types separetely. 
# Inputs: 
# 1. table-clean-filtered.qza (unrarified, decontam-filtered feature table that has been filtered for mitochondria and chloroplast)
# 2. taxonomy.qza
# 3. metadata.tsv (metadata without negative controls)
# Outputs: 
# 1. Volcano plot with W scores for each taxa that was found to be differentially abundant. There is a plot for mouth by area, cloaca by area, and feet by area. 



echo "Activating qiime environemnt"
source activate qiime2-2022.8

echo "starting sample filtering mouth" && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/dada2-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where "[type]='M'" --o-filtered-table ~/qiime/ANCOM/mouth-table.qza && \
echo "adding pseudocount to deal with zeros in table" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table.qza && \
echo "starting difference in mouth samples between two areas" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/ancom-area-mouth.qzv && \
echo "collapsing feature table to level 5" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/mouth-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 5 --o-collapsed-table ~/qiime/ANCOM/mouth-table-l5.qza && \
echo "add pseudocount" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table-l5.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table-l5.qza && \
echo "start composition differences at level 5" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table-l5.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l5-ancom-area-mouth.qzv && \
\
echo "starting sample filtering cloaca" && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/dada2-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where "[type]='C'" --o-filtered-table ~/qiime/ANCOM/cloaca-table.qza && \
echo "adding pseudocount to deal with zeros in table" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table.qza && \
echo "starting difference in mouth samples between two areas" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/ancom-area-cloaca.qzv && \
echo "collapsing feature table to level 5" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/cloaca-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 5 --o-collapsed-table ~/qiime/ANCOM/cloaca-table-l5.qza && \
echo "add pseudocount" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table-l5.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table-l5.qza && \
echo "start composition differences at level 5" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table-l5.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l5-ancom-area-cloaca.qzv && \
\
echo "starting sample filtering feet" && \
qiime feature-table filter-samples --i-table ~/qiime/denoising/dada2-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where "[type]='F'" --o-filtered-table ~/qiime/ANCOM/feet-table.qza && \
echo "adding pseudocount to deal with zeros in table" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/feet-table.qza --o-composition-table ~/qiime/ANCOM/comp-feet-table.qza && \
echo "starting difference in feet samples between two areas" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area  --o-visualization ~/qiime/ANCOM/ancom-area-feet.qzv && \
echo "collapsing feature table to level 5" && \
qiime taxa collapse --i-table ~/qiime/ANCOM/feet-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 5 --o-collapsed-table ~/qiime/ANCOM/feet-table-l5.qza && \
echo "add pseudocount" && \
qiime composition add-pseudocount --i-table ~/qiime/ANCOM/feet-table-l5.qza --o-composition-table ~/qiime/ANCOM/comp-feet-table-l5.qza && \
echo "start composition differences at level 5" && \
qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table-l5.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column Area --o-visualization ~/qiime/ANCOM/l5-ancom-area-feet.qzv && \
#\
#echo "start difference in feet samples between two body conditions"
#qiime feature-table filter-samples --i-table ~/qiime/ANCOM/bc-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --p-where "[type]='F'" --o-filtered-table ~/qiime/ANCOM/feet-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/ancom-bc-feet.qzv && \
#echo "collapsing feature table to level 6" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/feet-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/feet-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/feet-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-feet-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/l6-ancom-bc-feet.qzv && \
#\
#echo "start difference in cloaca samples between two body conditions"
#qiime feature-table filter-samples --i-table ~/qiime/denoising/dada2-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --p-where "[type]='C'" --o-filtered-table ~/qiime/ANCOM/cloaca-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/ancom-bc-cloaca.qzv && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/cloaca-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/cloaca-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/l6-ancom-bc-cloaca.qzv
#\
#echo "start difference in mouth samples between two body conditions"
#qiime feature-table filter-samples --i-table ~/qiime/denoising/dada2-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --p-where "[type]='F'" --o-filtered-table ~/qiime/ANCOM/feet-table.qza && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/ancom-bc-mouth.qzv && \
#echo "collapsing feature table to level 6" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/mouth-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/mouth-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata_rmbc.tsv --m-metadata-column condition_score --o-visualization ~/qiime/ANCOM/l6-ancom-bc-mouth.qzv && \
#\
#echo "start difference in feet samples between two sexes"
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/ancom-sex-feet.qzv && \
#echo "collapsing feature table to level 6" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/feet-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/feet-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/feet-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-feet-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-feet-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/l6-ancom-sex-feet.qzv && \
#\
#echo "start difference in cloaca samples between two body sexes"
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/ancom-sex-cloaca.qzv && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/cloaca-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/cloaca-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/cloaca-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-cloaca-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/l6-ancom-sex-cloaca.qzv
#\
#echo "start difference in mouth samples between two sexes"
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/ancom-sex-mouth.qzv && \
#echo "collapsing feature table to level 6" && \
#qiime taxa collapse --i-table ~/qiime/ANCOM/mouth-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-level 6 --o-collapsed-table ~/qiime/ANCOM/mouth-table-l6.qza && \
#echo "add pseudocount" && \
#qiime composition add-pseudocount --i-table ~/qiime/ANCOM/mouth-table-l6.qza --o-composition-table ~/qiime/ANCOM/comp-mouth-table-l6.qza && \
#echo "start composition differences at level 6" && \
#qiime composition ancom --i-table ~/qiime/ANCOM/comp-mouth-table-l6.qza --m-metadata-file ~/thesis/metadata/metadata_rmsex.tsv --m-metadata-column sex --o-visualization ~/qiime/ANCOM/l6-ancom-sex-mouth.qzv && \
date\

