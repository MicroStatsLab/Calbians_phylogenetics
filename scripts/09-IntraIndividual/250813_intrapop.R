library(tidyr)
library(readxl)
library(readr)
library(tidyverse)

TableS1 <- read_csv("Calb_manu_tables/TableS1.csv")
Intrapop_data <- meta_data %>%
  filter(Intrapop_ID != "NA")

TableS4 <- read_csv("Calb_manu_tables/TableS4.csv")

table(TableS4$Reference)
table(TableS4$Country)
table(TableS4$`Source of isolates`)
table(TableS4$`Sampling strategy`)

numSource <- table(TableS4$`Number of isolates`, TableS4$`Source of isolates`)
colSums(numSource)

nrow(Intrapop_data)
table(Intrapop_data$Intrapop_ID)
table(Intrapop_data$Intrapop_ID,Intrapop_data$Continent)


table(Intrapop_data$Intrapop_ID, Intrapop_data$Isolation_source)

# Create a frequency table of isolates per individual
isolate_counts <- table(Intrapop_data$Intrapop_ID)

# Total number of individuals
total_individuals <- length(isolate_counts)

# Mean number of isolates per individual
mean_isolates <- mean(isolate_counts)

# Median number of isolates per individual
median_isolates <- median(isolate_counts)

# Load dplyr for data manipulation

# Step 2: Get continent info (assuming each individual is from one continent)
individual_continent <- Intrapop_data %>%
  select(Intrapop_ID, Continent) %>%
  distinct()

# Step 3: Join isolate counts with continent info
isolate_by_continent <- isolate_counts %>%
  left_join(individual_continent, by = "Intrapop_ID")

# Step 4: Summarize: mean and median isolates per individual, by continent
summary_by_continent <- isolate_by_continent %>%
  group_by(Continent) %>%
  summarise(
    Total_Individuals = n(),
    Mean_Isolates = round(mean(Num_Isolates), 2),
    Median_Isolates = median(Num_Isolates),
    Max_Isolates = max(Num_Isolates),
    .groups = 'drop'
  )

# View results
print(summary_by_continent)


# Count number of unique clades per individual
clade_diversity <- Intrapop_data %>%
  group_by(Intrapop_ID) %>%
  summarise(n_clades = n_distinct(Proposed_final_clade)) %>%
  ungroup()

# Summary: how many individuals have multiple clades
table(clade_diversity$n_clades)

Intrapop_data_ddn <- split(Intrapop_data, Intrapop_data$Intrapop_ID)
t <- lapply(Intrapop_data_ddn, function(x) unique(x$Proposed_final_clade))
t <- lapply(Intrapop_data_ddn, function(x) unique(x$MAT_type))



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
Intrapop_data_with_frac_het <- Intrapop_data[!is.na(Intrapop_data$frac_HET), ]
# Count number of isolates per individual
isolate_counts <- table(Intrapop_data_with_frac_het$Intrapop_ID)

# Keep only individuals with more than 1 isolate
multi_isolate_ids <- names(isolate_counts[isolate_counts > 1])

# Filter the data
Intrapop_multi <- Intrapop_data_with_frac_het[Intrapop_data_with_frac_het$Intrapop_ID %in% multi_isolate_ids, ]

table(Intrapop_multi$Intrapop_ID)



####Visualize with boxplot
ggplot(Intrapop_multi, aes(x = Intrapop_ID, y = frac_HET)) +
  geom_boxplot(outlier.shape = NA, fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.4, color = "darkblue") +
  labs(
    title = "Fraction of Heterozygous Sites (frac_het) by Population",
    x = "Population",
    y = "Mean genome-wide heteroygosity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )



# # Now plot the filtered data
ggplot(Intrapop_multi, aes(x = Intrapop_ID, y = frac_HET)) +
  geom_violin(trim = TRUE, drop = FALSE) +  # Create the violin plot
  geom_sina() +  # Add individual data points
  labs(x = "Individual", y = "Mean genome-wide heteroygosity") +  # Axis labels
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size = 0.5, color = "red") +  # Add the mean line
  theme_minimal() +  # Clean theme
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  ylim(0, 0.012)



# Estimate the range of het perpopulation
fracHET_summary <- Intrapop_multi %>%
  group_by(Intrapop_ID) %>%
  summarise(
    N = n(),
    Min_fracHET = min(frac_HET, na.rm = TRUE),
    Max_fracHET = max(frac_HET, na.rm = TRUE),
    Mean_fracHET = mean(frac_HET, na.rm = TRUE),
    SD_fracHET = sd(frac_HET, na.rm = TRUE)
  ) %>%
  mutate(Range = paste0(round(Min_fracHET, 4), " – ", round(Max_fracHET, 4))) %>%
  select(Intrapop_ID, N, Range, Mean_fracHET, SD_fracHET)



view(fracHET_summary)

# Ensure Intrapop_ID matches the key in your metadata table
final_summary <- fracHET_summary %>%
  left_join(Intra_individual_meta, by = c("Intrapop_ID" = "Individual ID")) %>%
  select(Intrapop_ID, Sampling_strategy = `Sampling strategy`, 
         General_context = `General context`, 
         N, Range, Mean_fracHET, SD_fracHET)

view(final_summary)



###Now join het

biocontext_join <- Intrapop_multi %>%
  left_join(Intra_individual_meta, by = c("Intrapop_ID" = "Individual ID")) #%>%
select(Intrapop_ID, Sampling_strategy = `Sampling strategy`, 
       General_context = `General context`, 
       N, Range, Mean_fracHET, SD_fracHET)





# # Make sure Intrapop_ID is a factor ordered by sampling strategy
# biocontext_join <- biocontext_join %>%
#   arrange(`Sampling strategy`, Intrapop_ID) %>%
#   mutate(Intrapop_ID = factor(Intrapop_ID, levels = unique(Intrapop_ID)))

# Plot
ggplot(biocontext_join, aes(x = Intrapop_ID, y = frac_HET)) +
  geom_violin(trim = TRUE, fill = "lightblue", color = "black") +
  geom_sina(alpha = 0.6, size = 1, color = "darkblue") +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size = 0.5, color = "red") +
  labs(x = "Individual", y = "Mean genome-wide heterozygosity") +
  facet_wrap(~`Sampling strategy`, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, size = 9, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  ylim(0, 0.012)



# Make sure Intrapop_ID is a factor ordered by General context
biocontext_join <- biocontext_join %>%
  arrange(`General context`, Intrapop_ID) %>%
  mutate(Intrapop_ID = factor(Intrapop_ID, levels = unique(Intrapop_ID)))

# Plot
ggplot(biocontext_join, aes(x = Intrapop_ID, y = frac_HET)) +
  geom_violin(trim = TRUE, fill = "lightblue", color = "black") +
  geom_sina(alpha = 0.6, size = 1, color = "darkblue") +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size = 0.5, color = "red") +
  labs(x = "Individual", y = "Mean genome-wide heterozygosity") +
  facet_wrap(~`General context`, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, size = 9, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  ylim(0, 0.012)

