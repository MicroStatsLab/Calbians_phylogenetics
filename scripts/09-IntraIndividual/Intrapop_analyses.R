library(tidyverse)
library(here)
library(readxl)
library(tidyr)
library(car)
library(ggforce)
library(rstatix)
library(gridExtra)

# Load in data
TableS1 <- read_csv("Calb_manu_tables/TableS1.csv")
TableS1_intrapop <- TableS1 %>%
  filter(Intrapop_ID != "NA")

nrow(TableS1_intrapop)
table(TableS1_intrapop$Intrapop_ID)
table(TableS1_intrapop$Intrapop_ID,TableS1_intrapop$Continent)



table(TableS1_intrapop$Intrapop_ID, TableS1_intrapop$Isolation_source)

# Create a frequency table of isolates per individual
isolate_counts <- table(TableS1_intrapop$Intrapop_ID)

# Total number of individuals
total_individuals <- length(isolate_counts)

# Mean number of isolates per individual
mean_isolates <- mean(isolate_counts)

# Median number of isolates per individual
median_isolates <- median(isolate_counts)

individual_continent <- TableS1_intrapop %>%
  select(Intrapop_ID, Continent) %>%
  distinct()

# Count number of unique clades per individual
TableS1_intrapop_ddn <- split(TableS1_intrapop, TableS1_intrapop$Intrapop_ID)
t <- lapply(TableS1_intrapop_ddn, function(x) unique(x$Proposed_final_clade))

# multiple isolates
# IP17 - 15, 1 This study, Urine, Canada
# IP2 - 3, 1 Ford , oral, US
# IP22 - 4, 1, 15 Anderson, oral/Gastrointestinal, US
# IP28 - 15, 1 Anderson, oral, US
# IP34 - 3, 1 - Akhars, Oral, NA
# IP42 - 1, 2  - Akhars, Oral, NA
# IP43 1, 4  - Akhars, Oral, NA
# IP49 - 3, 1 - Akhars, Oral, NA
# IP5 - 2, 1 - Akhars, Oral, NA
# IP52 - 3, 2 - Akhars, Oral, NA
# IP59 - 1, 8 - Akhars, Oral, NA
# IP71 - 2, 8 - Akhars, Oral, NA
# IP76 - 4, 1 - Akhars, Oral, NA
# IP80 - 2, 4 - Akhars, Oral, NA
# IP91 = 3, 1 - Akhars, Oral, NA

#Extract intrapop isolates with frac_HET
TableS1_intrapop_with_frac_het <- TableS1_intrapop[!is.na(TableS1_intrapop$frac_HET), ]
# Count number of isolates per individual
isolate_counts_frac_het <- table(TableS1_intrapop_with_frac_het$Intrapop_ID)

# Keep only individuals with more than 1 isolate
multi_isolate_ids <- names(isolate_counts_frac_het[isolate_counts_frac_het > 1])

# Filter the data
intrapop_multi <- TableS1_intrapop_with_frac_het[TableS1_intrapop_with_frac_het$Intrapop_ID %in% multi_isolate_ids, ]

table(Intrapop_multi$Intrapop_ID)
