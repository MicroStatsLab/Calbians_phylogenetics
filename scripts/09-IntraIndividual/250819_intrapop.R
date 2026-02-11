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

TableS4 <- read_csv("Calb_manu_tables/TableS4.csv")
unique(TableS4$Reference)
#9 studies
table(TableS4$Sampling_strategy)
#Serial  Single Unknown 
#69      19       7 
table(TableS4$`General context`)
#Commensal   Disease 
#69        26 
table(TableS4$Continent)

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

# Estimate the range of het perpopulation
# fracHET_summary <- Intrapop_multi %>%
#   group_by(Intrapop_ID) %>%
#   summarise(
#     N = n(),
#     Min_fracHET = min(frac_HET, na.rm = TRUE),
#     Max_fracHET = max(frac_HET, na.rm = TRUE),
#     Mean_fracHET = mean(frac_HET, na.rm = TRUE),
#     SD_fracHET = sd(frac_HET, na.rm = TRUE)
#   ) %>%
#   mutate(Range = paste0(round(Min_fracHET, 4), " – ", round(Max_fracHET, 4))) %>%
#   select(Intrapop_ID, N, Range, Mean_fracHET, SD_fracHET)

# view(fracHET_summary)

# Ensure Intrapop_ID matches the key in your metadata table
# final_summary <- fracHET_summary %>%
#   left_join(TableS4, by = c("Intrapop_ID" = "Individual ID")) %>%
#   select(Intrapop_ID, Sampling_strategy = Sampling_strategy, 
#          General_context = `General context`, 
#          N, Range, Mean_fracHET, SD_fracHET)
# 
# view(final_summary)
# 


###Now join het
biocontext_join <- Intrapop_multi %>%
  left_join(TableS4, by = "Intrapop_ID") %>% 
  arrange(factor(Sampling_strategy))
#%>%
 # select(Intrapop_ID, Sampling_strategy = Sampling_strategy, 
#         General_context = `General context`, frac_HET
 #        N, Range, Mean_fracHET, SD_fracHET)

# # Make sure Intrapop_ID is a factor ordered by sampling strategy
# biocontext_join <- biocontext_join %>%
#   arrange(Sampling_strategy, Intrapop_ID) %>%
#   mutate(Intrapop_ID = factor(Intrapop_ID, levels = unique(Intrapop_ID)))

meanHet <- ggplot(biocontext_join, aes(x = Intrapop_ID, y = frac_HET, colour = Sampling_strategy)) +
    geom_point() +
    labs(x = "Individual", y = "Mean genome-wide heterozygosity") +
    #facet_wrap(~Sampling_strategy, scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(
      #axis.line = element_line(color = "black", linewidth = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 45, size = 9, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ylim(0, 0.012)
  
ggsave("Calb_manu_figures_tables/09-IntraIndividual/Intrapop_MeanHet.pdf", width = 7.5, height = 5)