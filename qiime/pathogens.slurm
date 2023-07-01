#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=alpha
#SBATCH --error=/home/aparsons/qiime/log/pathogen.log
#SBATCH --output=/home/aparsons/qiime/log/pathogen.log

#This slurm file makes a table of all gull pathogens from the paper "Pathogenic microorganisms associated with gulls and terns (Laridae)" (Hubalek 2021). 
#Input: feature table, taxonomy table, name of the pathogen as listed in taxonomy table, metadata file
#Output: .qzv file that tells you which samples have the pathogen and how many features were found in each sample


date && \
\
echo "Activating qiime environemnt" && \
\
source activate qiime2-2022.8 && \
echo "Start filtering" && \
\
echo"starting botulism" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_7; s__Clostridium_botulinum" --o-filtered-table ~/qiime/pathogens/table-botulism.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-botulism.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-botulism.qzv && \
\
echo"starting Eimeria praecox" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Eimeria_praecox" --o-filtered-table ~/qiime/pathogens/table-eimeria.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-eimeria.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-eimeria.qzv && \
\
echo"starting Mycobacterium" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales; f__Mycobacteriaceae; g__Mycobacterium; s__Mycobacterium_xenopi" --o-filtered-table ~/qiime/pathogens/table-mycobacterium.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-mycobacterium.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-mycobacterium.qzv && \
\
echo "starting Streptococcus pluranimalium" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_pluranimalium" --o-filtered-table ~/qiime/pathogens/table-Streptococcus-pluranimalium.qza && \ 
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Streptococcus-pluranimalium.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Streptococcus-pluranimalium.qzv && \
\
echo "starting Streptococcus suis" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_suis" --o-filtered-table ~/qiime/pathogens/table-Streptococcus-suis.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Streptococcus-suis.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Streptococcus-suis.qzv && \
\
echo "starting Streptococcus anginosus" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus_anginosus" --o-filtered-table ~/qiime/pathogens/table-Streptococcus-anginosus.qza
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Streptococcus-anginosus.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Streptococcus-anginosus.qzv && \
\
echo "starting Staphylococcus vitulinus" && \     
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Staphylococcales; f__Staphylococcaceae; g__Staphylococcus; s__Staphylococcus_vitulinus" --o-filtered-table ~/qiime/pathogens/table-Staphylococcus-vitulinus.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Staphylococcus-vitulinus.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Staphylococcus-vitulinus.qzv && \
\
echo "starting Staphylococcus equorum" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Staphylococcales; f__Staphylococcaceae; g__Staphylococcus; s__Staphylococcus_equorum" --o-filtered-table ~/qiime/pathogens/table-Staphylococcus-equorum.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Staphylococcus-equorum.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Staphylococcus-equorum.qzv && \
\
echo "starting Erysipelothrix_rhusiopathiae " && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Erysipelothrix; s__Erysipelothrix_rhusiopathiae" --o-filtered-table ~/qiime/pathogens/table-Erysipelothrix-rhusiopathiae.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Erysipelothrix-rhusiopathiae.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Erysipelothrix-rhusiopathiae.qzv && \
\
echo "starting Helicobacter_pullorum " && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Campilobacterota; c__Campylobacteria; o__Campylobacterales; f__Helicobacteraceae; g__Helicobacter; s__Helicobacter_pullorum" --o-filtered-table ~/qiime/pathogens/table-Helicobacter-pullorum.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Helicobacter-pullorum.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Helicobacter-pullorum.qzv && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Campylobacter-concisus.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Campylobacter-concisus.qzv && \
\
echo "starting Vibrio metschnikovii" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio; s__Vibrio_metschnikovii" --o-filtered-table ~/qiime/pathogens/table-Vibrio-metschnikovii.qza && \    
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Vibrio-metschnikovii.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Vibrio-metschnikovii.qzv && \
\
echo "starting Vibrio cholerae " && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio; s__Vibrio_cholerae" --o-filtered-table ~/qiime/pathogens/table-Vibrio-cholerae.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Vibrio-cholerae.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Vibrio-cholerae.qzv && \
\
echo "starting Riemerella anatipestifer" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Weeksellaceae; g__Riemerella; s__Riemerella_anatipestifer" --o-filtered-table ~/qiime/pathogens/table-Riemerella-anatipestifer.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Riemerella-anatipestifer.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Riemerella-anatipestifer.qzv && \
\
echo "starting Plesiomonas shigelloides" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Plesiomonas; s__Plesiomonas_shigelloides" --o-filtered-table ~/qiime/pathogens/table-Plesiomonas-shigelloides.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Plesiomonas-shigelloides.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Plesiomonas-shigelloides.qzv && \
\
echo "starting Mycoplasma zalophi" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Mycoplasmatales; f__Mycoplasmataceae; g__Mycoplasma; s__Mycoplasma_zalophi" --o-filtered-table ~/qiime/pathogens/table-Mycoplasma-zalophi.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Mycoplasma-zalophi.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Mycoplasma-zalophi.qzv && \
\
echo "starting Mycoplasma iowae" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Bacilli; o__Mycoplasmatales; f__Mycoplasmataceae; g__Mycoplasma; s__Mycoplasma_iowae" --o-filtered-table ~/qiime/pathogens/table-Mycoplasma-iowae.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-Mycoplasma-iowae.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-Mycoplasma-iowae.qzv && \
\
echo"starting clostridium perfrigens" && \
qiime taxa filter-table --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --p-include "p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_1; s__Clostridium_perfringens" --o-filtered-table ~/qiime/pathogens/table-clostridium.qza && \
\
qiime feature-table summarize --i-table ~/qiime/pathogens/table-clostridium.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/qiime/pathogens/table-clostridium.qzv && \
date\
