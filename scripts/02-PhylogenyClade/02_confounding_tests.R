library(readr)
library(rcompanion)
library(nnet)
library(brglm2)
library(here)
library(tidyverse)
library(vcd)
library(readxl)
###########################################
### read in and clean data
###########################################
TableS1_update <- read_excel(here("data_in", "Calb_manu_tables", "TableS1.xlsx"))
nrow(TableS1_update) #1178

#Select phylogenetically informative set
TableS1_update <- TableS1_update %>%
  filter(Intrapop_removed == 0)
nrow(TableS1_update) #945

TableS1_update <- TableS1_update %>%
  filter(Proposed_final_clade != 13)
nrow(TableS1_update) #908
 

#Check NA's 
sum(is.na(TableS1_update$Proposed_final_clade))
sum(is.na(TableS1_update$Isolation_source)) #17
sum(is.na(TableS1_update$Continent))


###Omit all NA's 
TableS1_update <- na.omit(TableS1_update[, c("Proposed_final_clade", "Isolation_source", "Continent")])
nrow(TableS1_update) #891

# Create contingency tables
# clade by source
table_clade_source <- table(TableS1_update$Proposed_final_clade,
                            TableS1_update$Isolation_source)
# clade by continent
table_clade_continent <- table(TableS1_update$Proposed_final_clade,
                               TableS1_update$Continent)
# source by continent
table_source_continent <- table(TableS1_update$Isolation_source, 
                                TableS1_update$Continent)

table_clade_source_continent <- table(TableS1_update$Proposed_final_clade, TableS1_update$Isolation_source, 
                                TableS1_update$Continent)

###########################################
### statistics
###########################################

# Cramér's V (see here https://www.ibm.com/docs/en/cognos-analytics/12.0.x?topic=terms-cramrs-v for interpretation)
cramer_clade_source     <- cramerV(table_clade_source) #0.2329 The fields are only weakly associated
cramer_clade_continent  <- cramerV(table_clade_continent) #0.3298 The fields are moderately associated.
cramer_source_continent <- cramerV(table_source_continent) #0.4889 The fields are moderately associated.

# Number of cells with counts less than 5
sum(table_clade_source < 5)
sum(table_clade_continent < 5)
sum(table_source_continent < 5)

# Number of cells with zero counts
sum(table_clade_source == 0)
sum(table_clade_continent == 0)
sum(table_source_continent == 0)

# Total number of cells in each table (for context)
length(table_clade_source)
length(table_clade_continent)
length(table_source_continent)

#Since standard chi-square assumptions are violated, run Monte Carlo simulation to approximate the p-value

# ---- Chi-square tests with simulated p-values (better for large tables) ----
chisq_clade_source     <- chisq.test(table_clade_source, simulate.p.value = TRUE, B = 100000)
chisq_clade_continent  <- chisq.test(table_clade_continent, simulate.p.value = TRUE, B = 100000)
chisq_source_continent <- chisq.test(table_source_continent, simulate.p.value = TRUE, B = 100000)

###########################################
### figures
###########################################

pdf("Figures/assoc_source_continent.pdf", width = 8, height = 6)
table_source_continent <- table_source_continent[, 
                                                 colnames(table_source_continent) != "NA"]

assoc(table_source_continent,shade = TRUE)
dev.off()

pdf("assoc_source_continent.pdf", width = 8, height = 6)
assoc(
  table_source_continent,
  main = "Association Between Isolation Source and Continent",
  shade = TRUE,
  labeling_args = list(
    rot_labels = c(45, 0),
    varnames = c("Source", "Continent"),
    gp_labels = gpar(fontsize = 9),
    gp_varnames = gpar(fontface = "bold", fontsize = 11)
  ),
  gp_shading = shading_Friendly,
  margins = c(5, 5, 4, 2)
)

dev.off()



######
TableS1_update <- read_excel(here("data_in", "Calb_manu_tables", "TableS1.xlsx"))


TableS1_update <- TableS1_update %>%
  filter(Intrapop_removed == 0)
nrow(TableS1_update)


TableS1_update <- TableS1_update %>%
  filter(Proposed_final_clade != 13)
nrow(TableS1_update)

# Filter only HSC isolates
hsc_df <- TableS1_update %>%
  filter(str_starts(Accession, "HSC"))  



###Omit all NA's 
hsc_df <- na.omit(hsc_df[, c("Proposed_final_clade", "Isolation_source", "Continent")])




# Create contingency table
contingency_table <- table(hsc_df$Proposed_final_clade, hsc_df$Isolation_source)

# -------------------------------------------
# 1. Check Expected Counts for Chi-squared
# -------------------------------------------
chi_check <- chisq.test(contingency_table)

# View expected frequencies
expected_counts <- chi_check$expected
summary(expected_counts < 5)  # Check for small expected values

# -------------------------------------------
# 2. Perform Chi-squared (with simulation if needed)
# -------------------------------------------
if(any(expected_counts < 5)) {
  message("Using Monte Carlo simulation due to low expected counts.")
  chi_test <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 100000)
} else {
  chi_test <- chi_check
}

# Print test result
print(chi_test)


