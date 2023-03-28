#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=taxonomy
#SBATCH --error=taxonomy.log
#SBATCH --output=taxonomy.log
	
date && \
\
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "starting taxonomy assignment" && \
qiime feature-classifier classify-sklearn --i-classifier TrainClassifier/classifier.qza --i-reads DADA2/rep-seqs-dada2.qza --verbose --o-classification taxonomy.qza && \
echo "finish taxonomy assignment" && \
\
echo "start visualization" && \
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv && \
echo "finished visualization" && \
\
echo "start barplot" && \
qiime taxa barplot --i-table DADA2/dada2-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.tsv --o-visualization TaxonomicAssignment/taxa-bar-plots.qzv && \
echo "finished barplot" && \
date\


