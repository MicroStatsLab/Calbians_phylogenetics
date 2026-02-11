
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)

X938_isolates_with_clade_13_MAT_meta_data <- read.csv("Calb_manu_tables/TableS1.csv")
filtered_data <- X938_isolates_with_clade_13_MAT_meta_data %>%
  filter(Intrapop_removed == 0, Proposed_final_clade != 13)#, Continent != "NA")

table(filtered_data$MAT_type,filtered_data$Continent)

filtered_data <- X938_isolates_with_clade_13_MAT_meta_data %>%
  filter(Intrapop_removed == 0, Proposed_final_clade != 13, Continent != "NA")

tab_continent <- table(filtered_data$MAT_type, filtered_data$Continent)
tab_continent_prop <- round(table(filtered_data$MAT_type, filtered_data$Continent)/nrow(filtered_data), 2)

# Fisher's Exact Test with simulated p-value
fisher.test(tab_continent, simulate.p.value = TRUE, B = 1e5)
#Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)

#data:  tab_continent
#p-value = 0.01505
#alternative hypothesis: two.sided
#prop.table(tab_continent, margin = 2)  # proportions within each continent

filtered_data_source <- X938_isolates_with_clade_13_MAT_meta_data %>%
  filter(Intrapop_removed == 0, Proposed_final_clade != 13, Isolation_source != "NA")

tab_source <- table(filtered_data_source$MAT_type, filtered_data_source$Isolation_source)
tab_continent_prop <- round(table(filtered_data$MAT_type, filtered_data$Continent)/nrow(filtered_data), 2)

fisher.test(tab_source, simulate.p.value = TRUE, B = 1e5)
# data:  tab_source
# p-value = 0.07548
# alternative hypothesis: two.sided

tab_clade <- table(filtered_data$MAT_type, filtered_data$Proposed_final_clade)
fisher.test(tab_clade, simulate.p.value = TRUE, B = 1e5)

#Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)
# 
# data:  tab_clade
# p-value = 0.01473
# alternative hypothesis: two.sided

# intrapop

intra <- X938_isolates_with_clade_13_MAT_meta_data[grep("IP", X938_isolates_with_clade_13_MAT_meta_data$Intrapop_ID), ]
# redo with full set
table(intra$MAT_type, intra$Intrapop_ID)
# diffent MAT, 4/95 populations:
#IP17, IP2, IP4, IP5

IP2 <- subset(X938_isolates_with_clade_13_MAT_meta_data, Intrapop_ID == "IP2")
#Ford, Oral, NA, due to chr5 aneuploidy

IP4 <- subset(X938_isolates_with_clade_13_MAT_meta_data, Intrapop_ID == "IP4")
#Ford, Oral, NA, due to chr5 aneuploidy

IP5 <- subset(X938_isolates_with_clade_13_MAT_meta_data, Intrapop_ID == "IP5")
#Ford, Oral, NA, due to chr5 aneuploidy

IP17 <- subset(X938_isolates_with_clade_13_MAT_meta_data, Intrapop_ID == "IP17")
#Ummmm? Why are so many HSCs intrapop? This is not right??
