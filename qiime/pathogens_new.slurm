#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=alpha
#SBATCH --error=/home/aparsons/qiime/log/pathogen_new.log
#SBATCH --output=/home/aparsons/qiime/log/pathogen_new.log

# This slurm file makes a table of all gull pathogens from the paper "Pathogenic microorganisms associated with gulls and terns (Laridae)" (Hubalek 2021).
# Input: feature table, taxonomy table, name of the pathogen as listed in taxonomy table, metadata file
# Output: .qzv file 


date && \
\
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
\
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_7; s__Clostridium_botulinum" , "p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Eimeria_praecox" , "p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales; f__Mycobacteriaceae; g__Mycobacterium; s__Mycobacterium_xenopi" , "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_pluranimalium" , "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_suis" , "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_anginosus" , "p__Firmicutes; c__Bacilli; o__Staphylococcales; f__Staphylococcaceae; g__Staphylococcus; s__Staphylococcus_vitulinus", "p__Firmicutes; c__Bacilli; o__Staphylococcales; f__Staphylococcaceae; g__Staphylococcus; s__Staphylococcus_equorum" , "p__Firmicutes; c__Bacilli; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Erysipelothrix; s__Erysipelothrix_rhusiopathiae" , "p__Campilobacterota; c__Campylobacteria; o__Campylobacterales; f__Helicobacteraceae; g__Helicobacter; s__Helicobacter_pullorum" , "p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio; s__Vibrio_metschnikovii" , "p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio; s__Vibrio_cholerae", "p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Weeksellaceae; g__Riemerella; s__Riemerella_anatipestifer", "p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Plesiomonas; s__Plesiomonas_shigelloides" , "p__Firmicutes; c__Bacilli; o__Mycoplasmatales; f__Mycoplasmataceae; g__Mycoplasma; s__Mycoplasma_zalophi" , "p__Firmicutes; c__Bacilli; o__Mycoplasmatales; f__Mycoplasmataceae; g__Mycoplasma; s__Mycoplasma_iowae" --o-filtered-table ~/qiime/pathogens/table-all-pathogens.qza && \
qiime feature-table summarize --i-table ~/qiime/pathogens/table-all-pathogens.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-all-pathogens.qzv && \
date \

