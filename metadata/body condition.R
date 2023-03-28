#library(tidyverse)

field_data = read.csv("~/thesis/metadata/field_data.csv")
head(field_data)

#subset morphometric measurements and mass
BC = select(field_data, c('TAR..mm.','CUL..mm.','SKULL..mm.','MN..mm.', 'MX..mm.','Mass..g.'))
head(BC)

#find which morphometric variable has the highest correlation with mass

#round(cor(BC),
#      digits = 2 # rounded to 2 decimals)
      
      
BC$Mass..g. = as.numeric(BC$Mass..g.)
BC$logmass = log(BC$Mass..g.)      
      
#Check for outliers in the plot

plot(BC$SKULL..mm., BC$Mass..g., main = "Scatter plot of Skull(mm) vs. Mass(g)", xlab = "Skull(mm)", ylab = "Mass(g)")

#find the mean skull length (105.49mm)

skull_mean = mean(BC$SKULL..mm.)
skull_mean

#find the standardized major axis slope (b = 3.167328)

library(smatr)

sma = sma(BC$Mass..g. ~ BC$SKULL..mm., BC, log = "xy", method= "SMA")
sma
plot(sma)

#Create an extra column in dataframe with body condition data

BC$body_condition = BC$Mass..g. * (skull_mean/BC$SKULL..mm.)^ 3.167328


#merge body condition column into field_data file
BC$Bird = c(1:65)
body_condition = subset(BC, select = c("Bird", "body_condition"))
merged_field_data = merge(field_data, body_condition, by = "Bird")

#Create 

write.table(merged_field_data, file = "~/thesis/metadata/metadata.csv", sep = ",", row.names = FALSE, na = "")
write.table(merged_field_data, file = "~/thesis/metadata/metadata.tsv", sep = "\t", row.names = FALSE, na = "")




