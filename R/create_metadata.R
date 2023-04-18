library(dplyr)
library(ggplot2)
library(tidyverse)
library(smatr)

field_data = read.csv("~/thesis/metadata/field_data.csv")
extraction_notes = read.csv("~/thesis/metadata/extraction_notes.csv")
argon_lab_data = read.csv("~/thesis/metadata/argon_lab_data.csv")

# Add column to record origin of sex data
field_data = rename(field_data, Sexing = Sex)
field_data$sex_origin = ifelse(field_data$sex == "", "", "Carly PCA")

# Convert sample ID to bird id
char_list = c("C", "F", "M")
extraction_notes$Sample_or_control = ifelse(grepl(paste(char_list, collapse = "|"), extraction_notes$SampleShort),
                                            "sample", "control")
extraction_notes$Bird = ifelse(extraction_notes$Sample_or_control == "sample",
                               substr(extraction_notes$SampleShort, 1, nchar(extraction_notes$SampleShort) - 1),
                               extraction_notes$SampleShort)
extraction_notes$Environmental_control = ifelse(extraction_notes$Sample_or_control == "control" & !grepl("X", extraction_notes$SampleShort),
                                                "true", "false")


# ------ SEX DATA -------

# Extract rows that have sex data
sex_data = subset(field_data, select = c("TAR", "CUL", "SKULL", "MN", "MX", "Mass", "sex", "Bird", "sex_origin"))
#sex_data = filter(sex_data, sex != "")
sex_data = filter(sex_data, Mass != "not taken")
sex_data = filter(sex_data, MX != "NA")
nrow(sex_data)

# Plot sample two dimensional graph
colors = ifelse(sex_data$sex == "F", "red",
                ifelse(sex_data$sex == "M", "blue", "black"))
plot(sex_data$TAR, sex_data$CUL, col=colors, pch=16)

# Do PCA on the sex data
dim_data = subset(sex_data, select = c("TAR", "CUL", "SKULL", "MN", "MX", "Mass"))
dim_data = apply(dim_data, 2, as.numeric)
pca_result <- prcomp(dim_data, scale=TRUE)
dim_data_pca = data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])
dim_data_pca$sex <- sex_data$sex
dim_data_pca$Bird <- sex_data$Bird
plot(dim_data_pca$PC1, dim_data_pca$PC2, col=colors, pch=16)

# Apply clustering to all data, based on PC1: f < 0 and m > 1
dim_data_pca$sex = ifelse(dim_data_pca$sex == "",
                          ifelse(dim_data_pca$PC1 < 0, "F",
                                 ifelse(dim_data_pca$PC1 > 1, "M", "")),
                          dim_data_pca$sex)

# Plot the results, show original data vs clustering data in different colors
dim_data_pca$sex_origin = sex_data$sex_origin
color_dict <- list(
  red = "Female",
  blue = "Male",
  pink = "Female predicted",
  lightblue = "Male predicted",
  black = "Unknown sex")

colors = ifelse(dim_data_pca$sex_origin == "Carly PCA",
                ifelse(dim_data_pca$sex == "F", "red", "blue"),
                ifelse(dim_data_pca$sex == "M", "lightblue",
                       ifelse(dim_data_pca$sex == "F", "pink", "black")))
plot(dim_data_pca$PC1, dim_data_pca$PC2, col=colors, pch=16, main="Sex clustering using PCA analysis")
color_list = c("red", "pink", "blue", "lightblue", "black")
color_meaning = sapply(color_list, function(x) color_dict[[x]])
legend("topright", legend = color_meaning, fill = color_list)

# Put sex data back into the original field data
sex_prediction = subset(dim_data_pca, select = c("Bird", "sex"))
field_data_without_sex = field_data[, !names(field_data) == "sex"]
field_data = merge(field_data_without_sex, sex_prediction, by = "Bird", all.x = TRUE)


# ------ BODY CONDITION INDEX -------

# subset morphometric measurements and mass
BC = select(field_data, c('TAR','CUL','SKULL','MN', 'MX','Mass', 'Bird'))

# Find which morphometric variable has the highest correlation with mass
BC = filter(BC, Mass != "not taken")
BC = filter(BC, MX != "NA")
BC = as.data.frame(sapply(BC, as.numeric))
round(cor(BC), digits = 2)
BC$logmass = log(BC$Mass)  

# Check for outliers in the plot
plot(BC$SKULL, BC$Mass, main = "Scatter plot of Skull(mm) vs. Mass(g)", xlab = "Skull(mm)", ylab = "Mass(g)")

# Find the mean skull length (105.49mm)
skull_mean = mean(BC$SKULL)

# Find the standardized major axis slope 
dim_BC = subset(BC, select = c("TAR", "CUL", "SKULL", "MN", "MX", "Mass"))
sma = sma(dim_BC$Mass ~ dim_BC$SKULL, dim_BC, log = "xy", method= "SMA")
plot(sma)

# Create an extra column with body condition data
BC$body_condition = BC$Mass * (skull_mean/BC$SKULL)^ 3.127186
field_data = merge(field_data, subset(BC, select = c("Bird", "body_condition")), by = "Bird", all.x = TRUE)


# Merge in field data for each experiment
all_data = merge(extraction_notes, field_data, by = "Bird", all.x = TRUE)
all_data = merge(all_data, argon_lab_data, by = "SampleShort", all.x = TRUE)

# Add sample id in the format that qiime expects
all_data = cbind("sample id" = all_data$SampleShort, all_data)
all_data$"sample id" = paste0("Sample-", all_data$"sample id")

# Filter out the environmental samples
all_data = subset(all_data, Environmental_control == "false")

#Filter out negative lab controls
all_data2 = subset(all_data, Sample_or_control == "sample")

# Write data back to csv and tsv file
write.table(all_data2, file = "~/thesis/metadata/metadata.csv", sep = ",", row.names = FALSE, na = "")
write.table(all_data2, file = "~/thesis/metadata/metadata.tsv", sep = "\t", row.names = FALSE, na = "")


# ------ PHYLOSEQ METADATA -------

# Also write metadata with special second row to make qza_to_phyloseq happy
# https://forum.qiime2.org/t/qiime2r-missing-sample/8681/24?page=2
phy_data = subset(all_data, select = c("sample id", "Bird", "barcodes", "TAR", "CUL", "SKULL", "MN", "MX", "Mass", "body_condition", "sex", "Area", "Colony", "eggs", "capture", "Sample_or_control", "type", "Concentration"))
new_row = ifelse(sapply(phy_data, is.numeric), "numerical", "categorical")
new_row[1] = "#q2:types"
phy_data = rbind(new_row, phy_data)
write.table(phy_data, file = "~/thesis/metadata/metadata_phyloseq.tsv", sep = "\t", row.names = FALSE, na = "", quote = FALSE)
