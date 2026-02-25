library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(here)

# Function to clean sample names
clean_sample_names <- function(name) {
  name <- gsub("^\\[\\d+\\]", "", name)  # Remove patterns like "[2]", "[3]", etc.
  name <- gsub("\\s+", "", name)         # Remove any spaces
  name <- gsub(":GT$", "", name)
  name <- gsub("_hapA", "", name) # Remove the ":GT" suffix
  return(name)
}

#############################################
# read in the raw genotype information and 
# convert from bases to binary 0/1 
# het/hom
#############################################
# Define file path and chunk size
# This comes from 
# bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' 250215_All_high_quality_samples.vcf.gz > 250215_All_high_quality_samples_genotypes.tsv

# Define bin size
bin_size <- 5000


# Read in the data
file_path1 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_1", "260113_genotypes_split_1.tsv.gz")
genotypes1 <- fread(file_path1, sep = "\t", header = TRUE)

# Convert genotypes to genetic states (1 = homozygous, 0 = heterozygous)
genotypes1 <- genotypes1 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

# Add a bin identifier
genotypes_with_bins_S1 <- genotypes1 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

# Count heterozygous SNPs per bin per sample
window_summary_S1 <- genotypes_with_bins_S1 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes1)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

# Calculate the start and end bin for each chromosome based on POS values
# Calculate the first and last bins for each chromosome
chromosome_bin_limits_S1 <- genotypes_with_bins_S1 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

# Find the global start and end bin across all chromosomes
chromosome_bin_limits_S1$CHROM <- gsub("_A", "", chromosome_bin_limits_S1$CHROM)
global_start_bin_S1 <- min(chromosome_bin_limits_S1$Start_Bin)
global_end_bin_S1 <- max(chromosome_bin_limits_S1$End_Bin)

final_table_S1 <- window_summary_S1 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S1, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isolSet1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 2
###########################

file_path2 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_1", "260113_genotypes_split_2.tsv.gz")
genotypes2 <- fread(file_path2, sep = "\t", header = TRUE)

genotypes2 <- genotypes2 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

genotypes_with_bins_S2 <- genotypes2 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S2 <- genotypes_with_bins_S2 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes2)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S2 <- genotypes_with_bins_S2 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S2$CHROM <- gsub("_A", "", chromosome_bin_limits_S2$CHROM)
global_start_bin_S2 <- min(chromosome_bin_limits_S2$Start_Bin)
global_end_bin_S2 <- max(chromosome_bin_limits_S2$End_Bin)

final_table_S2 <- window_summary_S2 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S2, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


###########################
# Isolate Set 3
###########################

file_path3 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_1", "260113_genotypes_split_3.tsv.gz")
genotypes3 <- fread(file_path3, sep = "\t", header = TRUE)

genotypes3 <- genotypes3 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

genotypes_with_bins_S3 <- genotypes3 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S3 <- genotypes_with_bins_S3 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes3)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S3 <- genotypes_with_bins_S3 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S3$CHROM <- gsub("_A", "", chromosome_bin_limits_S3$CHROM)
global_start_bin_S3 <- min(chromosome_bin_limits_S3$Start_Bin)
global_end_bin_S3 <- max(chromosome_bin_limits_S3$End_Bin)

final_table_S3 <- window_summary_S3 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S3, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 4
###########################

file_path_S4 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_2", "260113_genotypes_split_4.tsv.gz")
genotypes_S4 <- fread(file_path_S4, sep = "\t", header = TRUE)

genotypes_S4 <- genotypes_S4 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

genotypes_with_bins_S4 <- genotypes_S4 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S4 <- genotypes_with_bins_S4 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S4)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S4 <- genotypes_with_bins_S4 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S4$CHROM <- gsub("_A", "", chromosome_bin_limits_S4$CHROM)
global_start_bin_S4 <- min(chromosome_bin_limits_S4$Start_Bin)
global_end_bin_S4 <- max(chromosome_bin_limits_S4$End_Bin)

final_table_S4 <- window_summary_S4 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S4, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


###########################
# Isolate Set 5
###########################

file_path_S5 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_2", "260113_genotypes_split_5.tsv.gz")
genotypes_S5 <- fread(file_path_S5, sep = "\t", header = TRUE)

genotypes_S5 <- genotypes_S5 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

genotypes_with_bins_S5 <- genotypes_S5 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S5 <- genotypes_with_bins_S5 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S5)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S5 <- genotypes_with_bins_S5 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S5$CHROM <- gsub("_A", "", chromosome_bin_limits_S5$CHROM)
global_start_bin_S5 <- min(chromosome_bin_limits_S5$Start_Bin)
global_end_bin_S5 <- max(chromosome_bin_limits_S5$End_Bin)

final_table_S5 <- window_summary_S5 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S5, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 6
###########################

file_path_S6 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_2", "260113_genotypes_split_6.tsv.gz")
genotypes_S6 <- fread(file_path_S6, sep = "\t", header = TRUE)

genotypes_S6 <- genotypes_S6 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                        # Missing or malformed genotypes
  )))

genotypes_with_bins_S6 <- genotypes_S6 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S6 <- genotypes_with_bins_S6 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S6)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S6 <- genotypes_with_bins_S6 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S6$CHROM <- gsub("_A", "", chromosome_bin_limits_S6$CHROM)
global_start_bin_S6 <- min(chromosome_bin_limits_S6$Start_Bin)
global_end_bin_S6 <- max(chromosome_bin_limits_S6$End_Bin)

final_table_S6 <- window_summary_S6 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S6, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 7
###########################

file_path_S7 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_3", "260113_genotypes_split_7.tsv.gz")
genotypes_S7 <- fread(file_path_S7, sep = "\t", header = TRUE)

genotypes_S7 <- genotypes_S7 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                         # Missing or malformed genotypes
  )))

genotypes_with_bins_S7 <- genotypes_S7 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S7 <- genotypes_with_bins_S7 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S7)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S7 <- genotypes_with_bins_S7 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S7$CHROM <- gsub("_A", "", chromosome_bin_limits_S7$CHROM)
global_start_bin_S7 <- min(chromosome_bin_limits_S7$Start_Bin)
global_end_bin_S7 <- max(chromosome_bin_limits_S7$End_Bin)

final_table_S7 <- window_summary_S7 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S7, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 8
###########################

file_path_S8 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_3", "260113_genotypes_split_8.tsv.gz")
genotypes_S8 <- fread(file_path_S8, sep = "\t", header = TRUE)

genotypes_S8 <- genotypes_S8 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                         # Missing or malformed genotypes
  )))

genotypes_with_bins_S8 <- genotypes_S8 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S8 <- genotypes_with_bins_S8 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S8)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S8 <- genotypes_with_bins_S8 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S8$CHROM <- gsub("_A", "", chromosome_bin_limits_S8$CHROM)
global_start_bin_S8 <- min(chromosome_bin_limits_S8$Start_Bin)
global_end_bin_S8 <- max(chromosome_bin_limits_S8$End_Bin)

final_table_S8 <- window_summary_S8 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S8, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S8.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 9
###########################

file_path_S9 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_3", "260113_genotypes_split_9.tsv.gz")
genotypes_S9 <- fread(file_path_S9, sep = "\t", header = TRUE)

genotypes_S9 <- genotypes_S9 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                         # Missing or malformed genotypes
  )))

genotypes_with_bins_S9 <- genotypes_S9 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S9 <- genotypes_with_bins_S9 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S9)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S9 <- genotypes_with_bins_S9 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S9$CHROM <- gsub("_A", "", chromosome_bin_limits_S9$CHROM)
global_start_bin_S9 <- min(chromosome_bin_limits_S9$Start_Bin)
global_end_bin_S9 <- max(chromosome_bin_limits_S9$End_Bin)

final_table_S9 <- window_summary_S9 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S9, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S9.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

###########################
# Isolate Set 10
###########################

file_path_S10 <- here("data_out", "06-Heterozygosity", "fullVCF_zip_3", "260113_genotypes_split_10.tsv.gz")
genotypes_S10 <- fread(file_path_S10, sep = "\t", header = TRUE)

genotypes_S10 <- genotypes_S10 %>%
  mutate(across(3:ncol(.), ~ case_when(
    . == "0/0" ~ 1,                       # Homozygous reference
    grepl("^0/[1-9][0-9]*$", .) ~ 0,       # Heterozygous ref/alt (0/1, 0/20, 0/86)
    grepl("^[1-9][0-9]*/0$", .) ~ 0,       # Heterozygous alt/ref (1/0, 20/0)
    grepl("^[1-9][0-9]*/[1-9][0-9]*$", .) ~ 1, # Homozygous alternate (1/1, 39/39, 82/82)
    TRUE ~ NA_real_                         # Missing or malformed genotypes
  )))

genotypes_with_bins_S10 <- genotypes_S10 %>%
  mutate(Bin = (POS - 1) %/% bin_size + 1)

window_summary_S10 <- genotypes_with_bins_S10 %>%
  group_by(CHROM, Bin) %>%
  summarise(across(2:(length(genotypes_S10)-1), ~ sum(. == 0, na.rm = TRUE), .names = "het_count_{.col}"))

chromosome_bin_limits_S10 <- genotypes_with_bins_S10 %>%
  group_by(CHROM) %>%
  summarise(
    Start_Bin = (min(POS) - 1) %/% bin_size + 1,  # Calculate first bin
    End_Bin = (max(POS) - 1) %/% bin_size + 1     # Calculate last bin
  )

chromosome_bin_limits_S10$CHROM <- gsub("_A", "", chromosome_bin_limits_S10$CHROM)
global_start_bin_S10 <- min(chromosome_bin_limits_S10$Start_Bin)
global_end_bin_S10 <- max(chromosome_bin_limits_S10$End_Bin)

final_table_S10 <- window_summary_S10 %>%
  pivot_longer(
    cols = starts_with("het_count_"),
    names_to = "Sample",
    values_to = "Heterozygous_SNPs"
  ) %>%
  mutate(Sample = gsub("het_count_", "", Sample)) 

write.table(final_table_S10, here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S10.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


#########################
# combine the ten "final"
# tables
#########################

final_table_S1 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S1.tsv"), sep="\t", header=TRUE)
final_table_S2 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S2.tsv"), sep="\t", header=TRUE)
final_table_S3 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S3.tsv"), sep="\t", header=TRUE)
final_table_S4 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S4.tsv"), sep="\t", header=TRUE)
final_table_S5 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S5.tsv"), sep="\t", header=TRUE)
final_table_S6 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S6.tsv"), sep="\t", header=TRUE)
final_table_S7 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S7.tsv"), sep="\t", header=TRUE)
final_table_S8 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S8.tsv"), sep="\t", header=TRUE)
final_table_S9 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S9.tsv"), sep="\t", header=TRUE)
final_table_S10 <- read.table(here("data_out", "06-Heterozygosity", "GW_LOH_allSites", "260114_GW_LOH_isol_S10.tsv"), sep="\t", header=TRUE)

combined <- rbind(final_table_S1, final_table_S2, final_table_S3, final_table_S4, final_table_S5, final_table_S6, final_table_S7, final_table_S8, final_table_S9, final_table_S10)

# add the numeric phylogeny position information 
phyloPosition <- read.csv(here("data_out", "02-PhylogenyClade", 
                               "250813_908isolate_phylogeny_numericPosition.csv"))
combined <- left_join(combined, phyloPosition, by = "Sample")

# add required meta data
metadata <- read_csv("Calb_manu_tables/TableS1.csv")
metadata_sub <- metadata[, c("Accession","Isolate_Name", "Proposed_final_clade", "Continent", "Isolation_source")]

combined_meta <- merge(combined, metadata_sub, by = "Isolate_Name",by.x = "Sample", by.y = "Accession", all.x = TRUE)

# realized late in the game that one isolate with a high frequency of missing data was accidentally included in the LOH analysis. Remove it here. SRR11235412
combined_meta <- subset(combined_meta, Sample != "SRR11235412") #744
write.csv(combined_meta, "data_out/06-Heterozygosity/260114_744isolates_GW_LOH_table.csv", row.names = FALSE)
#250813_744isolates_GW_LOH_table.csv