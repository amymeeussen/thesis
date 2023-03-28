library(dplyr)
library(ggplot2)

field_data = read.csv("~/thesis/field_data.csv")
print("Header of field data:")
print(colnames(field_data))

# Add column to record orign of sex dat
field_data$sex_origin = ifelse(field_data$sex == "", "", "Carly PCA")


# Extract rows that have sex data
sex_data = subset(field_data, select = c("TAR..mm.", "CUL..mm.", "SKULL..mm.", "MN..mm.", "MX..mm.", "Mass..g.", "sex", "Bird..", "sex_origin"))
print("Header of sex data:")
colnames(sex_data)
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
dim_data_pca$Bird.. <- sex_data$Bird..
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
sex_prediction = subset(dim_data_pca, select = c("Bird..", "sex"))
fied_data_without_sex = field_data[, !names(field_data) == "sex"]
merged_fied_data = merge(fied_data_without_sex, sex_prediction, by = "Bird..")

# Write data back to csv file
write.table(merged_fied_data, file = "~/field_data_sexed.tsv", sep = "\t", row.names = TRUE)

