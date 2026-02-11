# --- Libraries ---
library(tidyverse)
library(Polychrome)
library(colorspace)
library(pophelper)
#--- File paths ---
qmatrix <- sprintf("data_out/03-Admixture/901_isolates/260109full_901_outfiles_K%d.qopt", c(6:20))

#qmatrix <- read.table("data_out/03-Admixture/full_isolate/260109full__outfiles_K6.qopt")

sample_ids <- read.table("data_out/03-Admixture/901_isolates/902_bam.txt")
qmatrix$ind <- sample_ids$V1
colnames(qmatrix)[1:6] <- paste0("K", 1:6)

q_long <- qmatrix %>%
  pivot_longer(
    cols = starts_with("K"),
    names_to = "Cluster",
    values_to = "Ancestry"
  )

GW_LOH_744isolates <- read.csv("data_out/06-Heterozygosity/legacy/250813_744isolates_GW_LOH_table.csv")
metadata <- GW_LOH_744isolates %>%
  distinct(Sample, Proposed_final_clade, PhyloPos) %>%  # keep unique rows
  arrange(PhyloPos)   

q_long <- left_join(q_long, metadata[, c("Sample", "Proposed_final_clade","PhyloPos")], by = c("ind" = "Sample"))


q_long$ind <- factor(q_long$ind, levels = metadata$PhyloPos)

colors20 <- hcl(
  h = seq(15, 375, length.out = 21)[-21],
  c = 50,   # moderate chroma → not deep
  l = 80    # high luminance → light
)


ggplot(q_long, aes(x = ind, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors20) +
  labs(x = "Sample (ordered by tree)", y = "Ancestry proportion", title = "ADMIXTURE K = 29") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )






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
