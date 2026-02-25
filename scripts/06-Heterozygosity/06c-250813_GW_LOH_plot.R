library(tidyverse)
library(here)
library(readxl)
library(tidyr)
library(car)
library(ggforce)
library(rstatix)
library(gridExtra)

####
# This runs the statistics script first
###

source(here("scripts", "06-Heterozygosity", "06b-heterozygosity_statistics.R"))


# This comes from 
# bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' 250215_All_738_high_quality_samples.vcf.gz > 250215_All_738_high_quality_samples_genotypes.tsv

TableS1 <- read_csv("Calb_manu_figures_tables/TableS1.csv")
GW_LOH_744isolates <- read.csv("data_out/06-Heterozygosity/250813_744isolates_GW_LOH_table.csv")

# select only the isolates that pass the HET filter
TableS1_LOH <- subset(TableS1, TableS1$Accession %in% unique(GW_LOH_744isolates$Sample))

GW_LOH_744isolates <- GW_LOH_744isolates %>%
  # Arrange by phylogenetic position
  arrange(PhyloPos) %>%
  # Reorder sample based on clade position
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Violin plot for genome-wide data by isolation source and clade

# Reorder by isolation_source  based on the frac_Het
# this is from the file "06c-heterozygosity_statistics.R                      
isoLevels <- cladeStatsSum$clade[order(cladeStatsSum$meanHet, decreasing = TRUE)]           

# # Now plot the filtered data
GW_LOH_isolation <- ggplot(TableS1_LOH, aes(x = Isolation_source, y = frac_HET)) +
  geom_violin(trim = TRUE, drop = FALSE) +  # Create the violin plot
  geom_sina() +  # Add individual data points
  labs(x = "Isolation source", y = "Mean genome-wide heteroygosity") +  # Axis labels
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, size = 0.5, color = "red") +  # Add the mean line
  theme_minimal() +  # Clean theme
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_discrete(limits = isoLevels) +
  ylim(0, 0.015)

# removes the 14 isolates that have NA for isolation source

# Save the plot as a PDF 
ggsave(here("Calb_manu_figures_tables","06-Heterozygosity", "250813_GW_LOH_track_Heterozysity_plot.pdf"), plot = GW_LOH_isolation, width = 7.5, height = 4, units = "in", dpi = 300)

# In keynote, add Tukey statistical letters from the HSD test in the heterozygosity_statistics script 

#################################################################################
# For the all isolate GW dataset, average by clade for comparison
#################################################################################
GW_LOH_744isolates_clade_ddn <- split(GW_LOH_744isolates, GW_LOH_744isolates$Proposed_final_clade)
GW_LOH_744isolates_cladeAVG_ddn <- lapply(GW_LOH_744isolates_clade_ddn, function(x) aggregate(x[c("Heterozygous_SNPs")], x[c("Bin", "CHROM", "Proposed_final_clade")], mean, na.rm=TRUE))

GW_LOH_744isolates_cladeAVG <- do.call(rbind, GW_LOH_744isolates_cladeAVG_ddn)


########################################################################
# Look at  heterozygosity across the genome in the outlier strains
########################################################################
lowStrains <- subset(TableS1_LOH, frac_HET < 0.00348) #3
highStrains <- subset(TableS1_LOH, frac_HET > 0.01) #7
outlierStrains <- rbind(lowStrains, highStrains)

outlierStrains_data <- subset(GW_LOH_744isolates, GW_LOH_744isolates$Sample %in% outlierStrains$Accession)

#there is a better way to do this but I don't know how
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR640897"] <- 0.00386127
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR6669928"] <- 0.00394297
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "HSC160"] <- 0.00347797
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "HSC157"] <- 0.00305421
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR1554310"] <- 0.00287436
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587743"] <- 0.01112309
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587701"] <- 0.01187267
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR1554304"] <- 0.01129727
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587853"] <- 0.01357767
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587820"] <- 0.01271985
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587818"] <- 0.01252991
outlierStrains_data$Het_ratio[outlierStrains_data$Sample == "SRR13587792"] <- 0.01390181

outlierStrains_data <- outlierStrains_data %>%
  # Arrange by clade position first
  arrange(desc(Het_ratio)) %>%
  # Reorder sample based on clade position
  mutate(Sample = factor(Sample, levels = unique(Sample)))


#adapt proposed final clade to get in the right order
# [1] SRR13587792 SRR13587853 SRR13587820 SRR13587818 SRR13587701 SRR1554304  SRR13587743
# [8] HSC160      HSC157      SRR1554310 
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587792"] <- 23
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587853"] <- 24
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587820"] <- 25
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587818"] <- 26
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587701"] <- 27
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR1554304"] <- 28
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR13587743"] <- 29
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "HSC160"] <- 31
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "HSC157"] <- 32
outlierStrains_data$fakeClade[outlierStrains_data$Sample == "SRR1554310"] <- 33

# create fake clade for plotting purposes in the average clade dataset
GW_LOH_744isolates_cladeAVG$Proposed_final_clade[GW_LOH_744isolates_cladeAVG$Proposed_final_clade == "S"] <- 22
GW_LOH_744isolates_cladeAVG$Proposed_final_clade <- as.numeric(GW_LOH_744isolates_cladeAVG$Proposed_final_clade)

GW_LOH_744isolates_cladeAVG$fakeClade <- ifelse(GW_LOH_744isolates_cladeAVG$Proposed_final_clade < 13, GW_LOH_744isolates_cladeAVG$Proposed_final_clade, GW_LOH_744isolates_cladeAVG$Proposed_final_clade-1)

# grab only the columns from outlier strains that match
outlierStrains_data_minimal <- outlierStrains_data[, names(GW_LOH_744isolates_cladeAVG)]

GW_LOH_744isolates_cladeAVGout <- rbind(GW_LOH_744isolates_cladeAVG, outlierStrains_data_minimal)


plot_LOH_GW_data <- ggplot(GW_LOH_744isolates_cladeAVGout, aes(x = Bin, y = fakeClade, fill = Heterozygous_SNPs)) +
  geom_tile() +  # Create a heatmap using tiles
  scale_y_reverse() +
  scale_fill_gradientn(
    colors = c("white", "red", "brown"),  # Define a custom color scale with brown for high values
    values = c(0, 0.5, 1),  # Define the relative positions of the colors
    limits = c(0, 150),  # Set the limits of the scale (up to 300)
    na.value = "white"  # Set the color for NA values
  ) +  
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(x = "Genome Bins", y = "", fill = "Heterozygous SNPs") +
  #labs(x = "Genome Bins", fill = "Heterozygous SNPs") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 1),  # Reduce y-axis text size to fit sample names
    panel.grid = element_blank()  # Remove grid lines for cleaner heatmap
  ) + 
  scale_x_continuous(
    expand = c(0,0)
  ) + 
  theme(panel.spacing = unit(0, "lines")) + 
  geom_vline(
    data = GW_LOH_744isolates_cladeAVG_ddn[[1]] %>% 
      group_by(CHROM) %>% 
      summarise(max_bin = max(Bin)), 
    aes(xintercept = max_bin), 
    color = "black", linewidth = 0.3)

# Save the plot as a PDF with the best dimensions
ggsave(here("Calb_manu_figures_tables","06-Heterozygosity", "250813_GW_LOH_cladeAVG_outliers.pdf"), plot = plot_LOH_GW_data, width = 13, height = 7, units = "in", dpi = 300)


#########################
# the outlier isolates 
############################
plot_LOH_GW_data_outliers <- ggplot(outlierStrains_data, aes(x = Bin, y = Sample, fill = Heterozygous_SNPs)) +
  geom_tile() +  # Create a heatmap using tiles
  scale_fill_gradientn(
    colors = c("white", "red", "brown"),  # Define a custom color scale with brown for high values
    values = c(0, 0.5, 1),  # Define the relative positions of the colors
    limits = c(0, 150),  # Set the limits of the scale (up to 300)
    na.value = "white"  # Set the color for NA values
  ) +  
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(x = "Genome Bins", y = "Samples", fill = "Heterozygous SNPs") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 3),  # Reduce y-axis text size to fit sample names
    panel.grid = element_blank()  # Remove grid lines for cleaner heatmap
  ) + 
  scale_x_continuous(
    expand = c(0,0)
  ) + 
  theme(panel.spacing = unit(0, "lines")) + 
  geom_vline(
    data = GW_LOH_744isolates_cladeAVG_ddn[[1]] %>% 
      group_by(CHROM) %>% 
      summarise(max_bin = max(Bin)), 
    aes(xintercept = max_bin), 
    color = "black", linewidth = 0.5)
#geom_hline(data = line_data, aes(yintercept = y_position), color = "black", linewidth = 0.3)

# Save the plot as a PDF with the best dimensions
ggsave(here("Calb_manu_figures_tables","06-Heterozygosity", "250813_GW_LOH_track_Heterozysity_plot_outliers.pdf"), plot = plot_LOH_GW_data_outliers, width = 8, height = 5, units = "in", dpi = 300)

##########
# Now plot all strains
##########

GW_LOH_744isolates <- GW_LOH_744isolates %>%
  # Arrange by phylogenetic position
  arrange(PhyloPos) %>%
  # Reorder sample based on clade position
  mutate(Sample = factor(Sample, levels = unique(Sample)))

###Create vector of the strain IDs that are the terminal strain for each clade - this was done visually
y_position <- c("SRR13587780", "SRR6669919", "SRR13587541", "SRR6669944", "SRR13587841","HSC097", "SRR13587646", "HSC116", "SRR13587634", "SRR13587524","SRR13587681", "SRR13587748", "SRR13587700", "SRR13587858", "SRR10209618", "HSC158", "SRR13741109", "SRR6669918", "SRR13587624", "SRR13587876")

# Convert the y_postion vector to a data frame
line_data <- data.frame(y_position = y_position)

plot_LOH_GW_data_full <- ggplot(GW_LOH_744isolates, aes(x = Bin, y = Sample, fill = Heterozygous_SNPs)) +
  geom_tile() +  # Create a heatmap using tiles
  scale_fill_gradientn(
    colors = c("white", "red", "brown"),  # Define a custom color scale with brown for high values
    values = c(0, 0.5, 1),  # Define the relative positions of the colors
    limits = c(0, 150),  # Set the limits of the scale (up to 300)
    na.value = "white"  # Set the color for NA values
  ) +  
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(x = "Genome Bins", y = "Samples", fill = "Heterozygous SNPs") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 1),  # Reduce y-axis text size to fit sample names
    panel.grid = element_blank()  # Remove grid lines for cleaner heatmap
  ) + 
  scale_x_continuous(
    expand = c(0,0)
  ) + 
  theme(panel.spacing = unit(0, "lines")) + 
  geom_vline(
    data = GW_LOH_744isolates_cladeAVG_ddn[[1]] %>% 
      group_by(CHROM) %>% 
      summarise(max_bin = max(Bin)), 
    aes(xintercept = max_bin), 
    color = "black", linewidth = 0.3) +
  geom_hline(data = line_data, aes(yintercept = y_position), color = "black", linewidth = 0.3)

# Save the plot as a PDF with the best dimensions
ggsave(here("Calb_manu_figures_tables","06-Heterozygosity", "250813_GW_LOH_track_Heterozysity_plot.pdf"), plot = plot_LOH_GW_data_full, width = 10, height = 8, units = "in", dpi = 300)

##################################################################
# Figure to show low and high positions across the genome
##################################################################

# This is the cen and telo positions, we don't want to include SNPs there in our analysis/figures
cenPos <- c(313, 1023, 1024, 1250, 1644, 1860, 2201, 2202, 2297, 2298, 2751, 2752)
TeloPos <- c(1, 2, 3, 638, 1085, 1445, 1766, 2005, 2212, 2402)


alldat <- data.frame(globalPos = seq_along(GW_LOH_744isolates_overall$Bin), (GW_LOH_744isolates_overall$Heterozygous_SNPs))
allDatCen <- subset(alldat, globalPos %in% c(cenPos-1, cenPos, cenPos+1))
names(allDatCen)[2] <- "Heterozygous_SNPs"
allDatTelo <- subset(alldat, globalPos %in% c(TeloPos-3, TeloPos-2, TeloPos-1, TeloPos, TeloPos+1, TeloPos+2))
names(allDatTelo)[2] <- "Heterozygous_SNPs"

# the low cutoff is log10(5) = less than 1 SNPs per 1000 bp (5 snps in a 5kb bin)

par(mfrow=c(2, 1), mar=c(0, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(seq_along(GW_LOH_744isolates_overall$Bin), log10(GW_LOH_744isolates_overall$Heterozygous_SNPs), xaxt="n", pch = ifelse(log10(GW_LOH_744isolates_overall$Heterozygous_SNPs) < log10(5), 19, 21), col = ifelse(log10(GW_LOH_744isolates_overall$Heterozygous_SNPs) < log10(5), "orange", "black"), yaxt="n", xlab="", ylab= "", cex=0.5, ylim=c(log10(0.01), log10(200)), xlim=c(100, 2770))
# none of the low points are on centromeres 
points(allDatCen$globalPos, log10(allDatCen$Heterozygous_SNPs), pch = ifelse(log10(allDatCen$Heterozygous_SNPs) < 0.6, 19, NA), col = ifelse(log10(allDatCen$Heterozygous_SNPs) < 0.6, "white", NA))
points(allDatTelo$globalPos, log10(allDatTelo$Heterozygous_SNPs), pch = ifelse(log10(allDatTelo$Heterozygous_SNPs) < 0.6, 19, NA), col = ifelse(log10(allDatTelo$Heterozygous_SNPs) < 0.6, "white", NA))
abline(v=638, lwd=0.5)
abline(v=638+447, lwd=0.5)
abline(v=638+447+360, lwd=0.5)
abline(v=638+447+360+321, lwd=0.5)
abline(v=638+447+360+321+239, lwd=0.5)
abline(v=638+447+360+321+239+207, lwd=0.5)
abline(v=638+447+360+321+239+207+190, lwd=0.5)
axis(1, at = c(638/2, 638+447/2, 638+447+360/2, 638+447+360+321/2, 638+447+360+321+239/2, 638+447+360+321+239+207/2, 638+447+360+321+239+207+190/2, 638+447+360+321+239+207+190+458/2), labels=FALSE)
axis(2, las=2, at= c(-1, 0, 1, 2), labels=c(0.1, 1, 10, 100))

# the high cutoff is 100 = more than than 20 SNPs per 1000 bp (100 snps in a 5kb bin)

plot(seq_along(GW_LOH_744isolates_overall$Bin), GW_LOH_744isolates_overall$Heterozygous_SNPs, xaxt="n", pch = ifelse(GW_LOH_744isolates_overall$Heterozygous_SNPs > 100, 19, 21), cex= 0.5, col = ifelse(GW_LOH_744isolates_overall$Heterozygous_SNPs > 100, "purple", "black"), yaxt="n", xlab="", ylab= "", ylim=c(0, 200), xaxs="i", xlim=c(-10, 2870))
points(allDatCen$globalPos, allDatCen$Heterozygous_SNPs, pch = ifelse(allDatCen$Heterozygous_SNPs > 100, 19, NA), col = ifelse(allDatCen$Heterozygous_SNPs > 100, "white", NA))
points(allDatTelo$globalPos, allDatTelo$Heterozygous_SNPs, pch = ifelse(allDatTelo$Heterozygous_SNPs > 100, 19, NA), col = ifelse(allDatTelo$Heterozygous_SNPs > 100, "white", NA))
#arrows(cenPos, rep(200, length(cenPos)), cenPos, rep(0, length(cenPos)),pch=4, cex=0.7, col="grey", length=0)
abline(v=638, lwd=0.5)
abline(v=638+447, lwd=0.5)
abline(v=638+447+360, lwd=0.5)
abline(v=638+447+360+321, lwd=0.5)
abline(v=638+447+360+321+239, lwd=0.5)
abline(v=638+447+360+321+239+207, lwd=0.5)
abline(v=638+447+360+321+239+207+190, lwd=0.5)
#abline(v=638+447+360+321+239+207+190+458, lwd=0.5)
axis(1, at = c(638/2, 638+447/2, 638+447+360/2, 638+447+360+321/2, 638+447+360+321+239/2, 638+447+360+321+239+207/2, 638+447+360+321+239+207+190/2, 638+447+360+321+239+207+190+458/2), labels=c(paste0("chr", 1:7), "chrR"))
axis(2, las=2)
mtext("Average number of SNPs (/5 kb bin)", side=2, 
      outer=TRUE, line=2)

