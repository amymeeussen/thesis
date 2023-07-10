# This RScript reads the metadata file and generates a histogram plot of the body conditions
# in the San Francisco and Mono lake locations. Also tests significance of difference between
# body condition in the two locations.

# Input: ~/thesis/metadata/metadata.csv file on disk
# Output: histogram plot


library(dplyr)
metadata = read.csv("~/thesis/metadata/metadata.csv")

blue = rgb(83/255,28/255,136/255)
green = rgb(148/255,209/255,35/255)

ml_body_condition = metadata %>%
  filter(Area == "ML") %>%
  select(body_condition)
sf_body_condition = metadata %>%
  filter(Area == "SF") %>%
  select(body_condition)

# Test significance
result = t.test(na.omit(sf_body_condition), na.omit(ml_body_condition), var.equal = FALSE)


# Create histograms
ml_hist = hist(ml_body_condition$body_condition, breaks = seq(350, 850, by = 50))
sf_hist = hist(sf_body_condition$body_condition, breaks = seq(350, 850, by = 50))

ml_df <- data.frame(
  breaks = ml_hist$breaks[-length(ml_hist$breaks)],
  counts = ml_hist$counts
)
sf_df <- data.frame(
  breaks = sf_hist$breaks[-length(sf_hist$breaks)],
  counts = sf_hist$counts
)
df = merge(ml_df, sf_df, by="breaks")
barplot(
  rbind(df$counts.x, df$counts.y),
  beside = TRUE,
  names.arg = df$breaks,
  legend.text = c("Mono Lake", "San Francisco"),
  col = c(blue, green),
  xlab = "Body Condition(scaled mass index)",
  ylab = "Number of birds", 
  cex.axis = 1.5,
  cex.names = 1.5, 
  cex.lab = 1.5) 
