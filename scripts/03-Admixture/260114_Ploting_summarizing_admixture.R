# --- Libraries ---
library(tidyverse)
library(Polychrome)
library(colorspace)
library(pophelper)
#--- File paths ---
qmatrix <- sprintf("data_out/03-Admixture/901_isolates/260109full_901_outfiles_K%d.qopt", c(6:20))



#Calculating the number of ancestries per isolate
# K values we are interested in
Ks <- c(12, 16, 20)

# Set a threshold for “ancestry present”
threshold <- 0.05

summarize_ancestries <- function(K, sample_file = "data_out/03-Admixture/901_isolates/902_bam.txt") {
  
  # Read Q matrix
  qfile <- sprintf("data_out/03-Admixture/901_isolates/260109full_901_outfiles_K%d.qopt", K)
  qmat <- read.table(qfile)
  
  # Add sample IDs
  sample_ids <- read.table(sample_file)
  qmat$Sample <- sample_ids$V1
  
  # Count number of ancestries above threshold for each sample
  ancestry_count <- qmat %>%
    rowwise() %>%
    mutate(n_ancestries = sum(c_across(1:K) > threshold)) %>%
    select(Sample, n_ancestries)
  
  # Add K column
  ancestry_count$K <- K
  
  return(ancestry_count)
}

result <- bind_rows(lapply(Ks, summarize_ancestries))

library(tidyr)

result_wide <- result %>%
  pivot_wider(names_from = K, values_from = n_ancestries, names_prefix = "K")

result_wide
View(result_wide)


metadata_1 <- read_csv("Calb_manu_tables/TableS1.csv")

result_wide <- result_wide %>%
  left_join(metadata_1 %>% select(Accession, Proposed_final_clade,Isolate_Name),
            by = c("Sample" = "Accession"))


View(result_wide)

write.table(result_wide_1, "data_out/03-Admixture/901_isolates/260114_ancestries_prop.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#############
AncProp <- read.csv("data_out/03-Admixture/901_isolates/260114_ancestries_prop.tsv", sep = "\t")
TableS1_admix <- read.csv("data_out/03-Admixture/901_isolates/TableS1_admix.csv")
ad <- subset(TableS1_admix, K20 > 0)

barplot(table(ad$K20, ad$Isolation_source))
ad_high <- subset(ad, K20 > 6)
table(ad_high$Isolation_source)
plot(ad$K20, ad$frac_HET)

highHet <- subset(ad, frac_HET > 0.01)
