
field_data = read.csv("~/thesis/metadata/field_data.csv", na.strings = c("", "NA"))
extraction_notes = read.csv("~/thesis/metadata/extraction_notes.csv", na.strings = c("", "NA"))
argon_lab_data = read.csv("~/thesis/metadata/argon_lab_data.csv", na.strings = c("", "NA"))

# Convert sample ID to bird id
extraction_notes$Bird = ifelse(grepl("\\.", extraction_notes$Sample.ID),
                               extraction_notes$Sample.ID,
                               substr(extraction_notes$Sample.ID, 1, nchar(extraction_notes$Sample.ID) - 1))

# Merge in field data for each experiment
all_data = merge(extraction_notes, field_data, by = "Bird", all.x = TRUE)
all_data = merge(all_data, argon_lab_data, by = "Sample.ID", all.x = TRUE)

# Write the end result to file
write.table(all_data, file = "~/thesis/metadata/metadata.csv", sep = ",", row.names = FALSE, na = "")
