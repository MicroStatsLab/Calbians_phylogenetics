library(dplyr)
library(tidyr)
library(purrr)
library(readr)
K <- 6:20

q_long_all <- map_dfr(K, function(K) {
  
  qmatrix <- read.table(
    paste0("data_out/03-Admixture/901_isolates/260109full_901_outfiles_K", K, ".qopt")
  )
  
  sample_ids <- read.table("data_out/03-Admixture/901_isolates/902_bam.txt")
  
  qmatrix$ind <- sample_ids$V1
  colnames(qmatrix)[1:K] <- paste0("K", 1:K)
  
  qmatrix %>%
    pivot_longer(
      cols = starts_with("K"),
      names_to = "Cluster",
      values_to = "Ancestry"
    ) %>%
    mutate(
      K = K   # 🔑 THIS is what was missing
    )
})
sorted_df <- read_csv("data_out/02-PhylogenyClade/250813_908isolate_phylogeny_numericPosition.csv")
TableS1 <- read_csv("Calb_manu_figures_tables/TableS1.csv")


sorted_df <- sorted_df %>%
  left_join(
    TableS1 %>% select(Accession, Proposed_final_clade),
    by = c("Sample" = "Accession")
  )

q_long_all$ind <- factor(q_long_all$ind, levels = sorted_df$Sample)

clade_boundaries <- sorted_df %>%
  mutate(pos = row_number()) %>%
  group_by(Proposed_final_clade) %>%
  summarise(end_pos = max(pos) + 0.5, .groups = "drop")


# cols <- c(
#   "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
#   "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5",
#   "#FFED6F", "#A6CEE3", "#FDBF6F", "#CAB2D6", "#E6F5C9",
#   "#FDE0DD", "#D0E1F2", "#FEE8C8", "#E5C1CD", "#CDE5D8"
# )

cols <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#F46D43", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3DE69", "#FDBF6F", "#80B1D3", "#FB8072", "#BEBADA"
)

q_long_all <- q_long_all %>%
  mutate(
    ind = factor(ind, levels = sorted_df$Sample),
    K   = factor(K, levels = 6:20, ordered = TRUE)
  )


ggplot(q_long_all, aes(x = ind, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  geom_vline(
    data = clade_boundaries,
    aes(xintercept = end_pos),
    colour = "black",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ K, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols) +
  labs(
    x = "Sample (ordered by tree)",
    y = "Ancestry proportion"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",     # ⬅️ exclude legend
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.01, "lines"),
    strip.text = element_text(face = "bold",size=7)
  )



####
sample_pos <- sorted_df %>%
  mutate(pos = row_number()) %>%
  select(Sample, pos)

q_long_all <- q_long_all %>%
  left_join(sample_pos, by = c("ind" = "Sample")) %>%
  mutate(
    K = factor(K, levels = 6:20, ordered = TRUE)
  )

clade_boundaries <- sorted_df %>%
  mutate(pos = row_number()) %>%
  group_by(Proposed_final_clade) %>%
  summarise(end_pos = max(pos), .groups = "drop")


ggplot(q_long_all, aes(x = pos, y = Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  geom_vline(
    data = clade_boundaries,
    aes(xintercept = end_pos),
    colour = "black",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ K, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Sample (ordered by tree)",
    y = "Ancestry proportion"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.01, "lines"),
    strip.text = element_text(face = "bold", size = 7)
  )


#####Strangely i have to run this twice

library(dplyr)
library(ggplot2)
library(grid)

# -------------------------
# 1️⃣ Keep only non-empty samples
# -------------------------
nonempty_samples <- q_long_all %>%
  group_by(ind) %>%
  summarise(total_ancestry = sum(Ancestry), .groups = "drop") %>%
  filter(total_ancestry > 0) %>%
  pull(ind)

q_long_all <- q_long_all %>%
  filter(ind %in% nonempty_samples)

# -------------------------
# 2️⃣ Assign positions based only on non-empty samples
# -------------------------
sample_pos <- sorted_df %>%
  filter(Sample %in% nonempty_samples) %>%
  mutate(pos = row_number()) %>%
  select(Sample, pos, Proposed_final_clade)

q_long_all <- q_long_all %>%
  left_join(sample_pos, by = c("ind" = "Sample")) %>%
  mutate(K = factor(K, levels = 6:20, ordered = TRUE))

# Make 'ind' a factor sorted by pos
q_long_all <- q_long_all %>%
  mutate(ind = factor(ind, levels = sample_pos$Sample))

# -------------------------
# 3️⃣ Clade boundaries
# -------------------------
clade_boundaries <- sample_pos %>%
  group_by(Proposed_final_clade) %>%
  summarise(end_pos = max(pos), .groups = "drop")

# -------------------------
# 4️⃣ Plot
# -------------------------
p <- ggplot(q_long_all, aes(x = pos, y = Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  geom_vline(data = clade_boundaries,
             aes(xintercept = end_pos),
             colour = "black",
             linewidth = 0.3) +
  facet_wrap(~ K, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  labs(x = "Sample (ordered by tree)", y = "Ancestry proportion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.01, "lines"),
    strip.text = element_text(face = "bold", size = 7)
  )


# Save as PDF (vector graphics, maximum quality)
ggsave(
  filename = "ancestry_plot.pdf",  # output file
  plot = p,                        # your ggplot object
  width = 14,                      # width in inches
  height = 13,                      # height in inches
  units = "in"                      # units for width/height
)

