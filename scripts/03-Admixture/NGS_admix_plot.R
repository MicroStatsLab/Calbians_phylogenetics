# --- Libraries ---
library(tidyverse)
library(Polychrome)
library(colorspace)
library(pophelper)

#' #Analyses run on fir 
#' #Generate the beagle file
#' #cd /home/abdul/scratch/alinmet_stats/BAM_files/passed_samples_final/NGSadmix
#' # Extract marker, alleles, and per-sample PLs from a biallelic VCF,
#' # then convert PL (phred-scaled likelihoods) to GL (genotype likelihoods)
#' # in BEAGLE format.
#' 
#' #bcftools query -f '%CHROM:%POS\t%REF\t%ALT[\t%PL]\n' vcf_biallelic.vcf.gz | \
#' #perl -F'\t' -ane '
#'     # Make output with BEAGLE-required fields:
#'     # marker ID, reference allele, alternate allele
#'     @out = ($F[0], $F[1], $F[2]);
#' 
#'     # Loop over samples (each PL field)
#'     for ($i = 3; $i < @F; $i++) {
#' 
#'         # Split comma-separated PLs
#'         @pl = split(",", $F[$i]);
#' 
#'         # Ensure exactly 3 PLs (biallelic: 0/0, 0/1, 1/1)
#'         if (@pl != 3) { @pl = (".", ".", "."); }
#' 
#'         # Convert PL → GL
#'         # GL = 10^(-PL/10); missing PLs become 0
#'         @gl = map { ($_ eq "." ? 0 : 10 ** (-$_/10)) } @pl;
#' 
#'         # If all GLs are zero (completely missing information),
#'         # replace with uniform probabilities
#'         $sum = 0; $sum += $_ for @gl;
#'         if ($sum == 0) { @gl = (1/3, 1/3, 1/3); }
#' 
#'         # Append genotype likelihoods to output
#'         push @out, @gl;
#'     }
#' 
#'     # Print BEAGLE-formatted line
#'     print join("\t", @out), "\n";
#' ' > output.beagle
#' 
#' 
#' # Quality checks: ensure the number of columns matches the header
#' awk -F'\t' '{print NF; exit}' output.beagle
#' awk -F'\t' '{print NF; exit}' header.txt
#' 
#' 
#' # Combine header and data into final BEAGLE file
#' cat header.txt output.beagle > 251125_final.beagle
#' #cat header.txt output.beagle > 251125_final.beagle
#' # Run this on the HPC
#' # module load angsd samtools bcftools htslib 
#' # # Run NGSadmix for K = 1 to 40 using GNU parallel
#' # seq 1 40 | parallel -j 10  \
#' # NGSadmix -likes 251125_final.beagle -K {} -P 32 -o outfiles_K{} -minMaf 0.05
#' 

 #--- File paths ---
qfiles <- sprintf("/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture/901_isolates/260109full_901_outfiles_K%d.qopt", c(1:15))

qfiles <- sprintf("/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture/full_isolate/260109full__outfiles_K%d.qopt", c(8:15))


###Admixture plot with K=10
qmatrix <- read.table("data_out/03-Admixture/full_isolate/260109full__outfiles_K6.qopt")

sample_ids <- read.table("data_out/03-Admixture/full_isolate/pop.list")
qmatrix$ind <- sample_ids$V1
colnames(qmatrix)[1:6] <- paste0("K", 1:6)


q_long <- qmatrix %>%
  pivot_longer(
    cols = starts_with("K"),
    names_to = "Cluster",
    values_to = "Ancestry"
  )


unique_samples <- GW_LOH_744isolates[!duplicated(GW_LOH_744isolates$Sample),
                                     c("Sample", "PhyloPos","Proposed_final_clade")]

sorted_df <- unique_samples %>%
  arrange(PhyloPos)




q_long$ind <- factor(q_long$ind, levels = sorted_df$Sample)



ggplot(q_long, aes(x = ind, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  labs(x = "Sample (ordered by tree)", y = "Ancestry proportion", title = "ADMIXTURE K = 2") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



library(dplyr)

clade_boundaries <- sorted_df %>%
  mutate(pos = row_number()) %>%
  group_by(Proposed_final_clade) %>%
  summarise(end_pos = max(pos), .groups = "drop")



ggplot(q_long, aes(x = ind, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  geom_vline(
    data = clade_boundaries,
    aes(xintercept = end_pos + 0.5),
    colour = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Sample (ordered by tree)",
    y = "Ancestry proportion",
    title = "ADMIXTURE K = 6"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Load metadata and sample order only once
sample_ids <- read.table("/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture/full_isolate/pop.list", stringsAsFactors = FALSE)
#metadata <- read.table("data_out/03_population_structure/sample_metadata.txt", header = TRUE, sep = "\t")  # adjust path and separator
phylo_order <- read.table("data_out/04_admixture/phylogeny_position.txt", stringsAsFactors = FALSE)




# --- Read Q-matrices into a qlist ---
qlist <- readQ(qfiles)

cols <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#FDBF6F", "#CAB2D6"
)

cols20 <- hcl.colors(20, "Set3")

# --- Align clusters across Ks to fix label switching ---
qlist_aligned <- alignK(qlist, type = "across")
plotQ(qlist_aligned,imgoutput="join",clustercol = cols ,exportpath="/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture/full_isolate")
