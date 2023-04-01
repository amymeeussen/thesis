library(dplyr)
library(ggplot2)
library(tidyverse)
library(smatr)

field_data = read.csv("~/thesis/metadata/field_data.csv")
extraction_notes = read.csv("~/thesis/metadata/extraction_notes.csv")
argon_lab_data = read.csv("~/thesis/metadata/argon_lab_data.csv")

# Add column to record origin of sex data
field_data$sex_origin = ifelse(field_data$sex == "", "", "Carly PCA")

# Convert sample ID to bird id
extraction_notes$Bird = ifelse(grepl("\\.", extraction_notes$Sample.ID),
                               extraction_notes$Sample.ID,
                               substr(extraction_notes$Sample.ID, 1, nchar(extraction_notes$Sample.ID) - 1))

# ------ SEX DATA -------

# Extract rows that have sex data
sex_data = subset(field_data, select = c("TAR..mm.", "CUL..mm.", "SKULL..mm.", "MN..mm.", "MX..mm.", "Mass..g.", "sex", "Bird", "sex_origin"))
#sex_data = filter(sex_data, sex != "")
sex_data = filter(sex_data, Mass..g. != "not taken")
sex_data = filter(sex_data, MX..mm. != "NA")
nrow(sex_data)

# Plot sample two dimensional graph
colors = ifelse(sex_data$sex == "F", "red",
                ifelse(sex_data$sex == "M", "blue", "black"))
plot(sex_data$TAR..mm., sex_data$CUL..mm., col=colors, pch=16)

# Do PCA on the sex data
dim_data = subset(sex_data, select = c("TAR..mm.", "CUL..mm.", "SKULL..mm.", "MN..mm.", "MX..mm.", "Mass..g."))
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
BC = select(field_data, c('TAR..mm.','CUL..mm.','SKULL..mm.','MN..mm.', 'MX..mm.','Mass..g.', 'Bird'))

# Find which morphometric variable has the highest correlation with mass
BC = filter(BC, Mass..g. != "not taken")
BC = filter(BC, MX..mm. != "NA")
BC <- as.data.frame(sapply(BC, as.numeric))
round(cor(BC), digits = 2)
BC$logmass = log(BC$Mass..g)  

# Check for outliers in the plot
plot(BC$SKULL..mm., BC$Mass..g., main = "Scatter plot of Skull(mm) vs. Mass(g)", xlab = "Skull(mm)", ylab = "Mass(g)")

# Find the mean skull length (105.49mm)
skull_mean = mean(BC$SKULL..mm.)

# Find the standardized major axis slope 
dim_BC = subset(BC, select = c("TAR..mm.", "CUL..mm.", "SKULL..mm.", "MN..mm.", "MX..mm.", "Mass..g."))
sma = sma(dim_BC$Mass..g. ~ dim_BC$SKULL..mm., dim_BC, log = "xy", method= "SMA")
plot(sma)

# Create an extra column with body condition data
BC$body_condition = BC$Mass..g. * (skull_mean/BC$SKULL..mm.)^ 3.127186
field_data = merge(field_data, BC, by = "Bird", all.x = TRUE)


# Merge in field data for each experiment
all_data = merge(extraction_notes, field_data, by = "Bird", all.x = TRUE)
all_data = merge(all_data, argon_lab_data, by = "Sample.ID", all.x = TRUE)

# Add sample id in the format that qiime expects
all_data = cbind("sample id" = all_data$Sample.ID, all_data)
all_data$"sample id" = paste0("Sample-", all_data$"sample id")

# Write data back to csv and tsv file
write.table(all_data, file = "~/thesis/metadata/metadata.csv", sep = ",", row.names = FALSE, na = "")
write.table(all_data, file = "~/thesis/metadata/metadata.tsv", sep = "\t", row.names = FALSE, na = "")
