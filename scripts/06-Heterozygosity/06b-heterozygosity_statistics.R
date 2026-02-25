library(tidyverse)
library(here)
library(readxl)
library(tidyr)
library(car)
library(ggforce)
library(rstatix)
library(gridExtra)
library(agricolae)

####################################
# Look at average heterozygosity across the genome
####################################

# The number of HOM and HET sites and chromosome stats from each isolate were calculated in plink 
# 250216_merged_WGSandchr_plink_data.tsv
# and added to Table S1.

TableS1 <- read_csv("Calb_manu_figures_tables/TableS1.csv")

GW_LOH_744isolates <- read.csv("data_out/06-Heterozygosity/GW_LOH_tables/260114_744isolates_GW_LOH_table.csv")

# select only the isolates that pass the HET filter
TableS1_LOH <- subset(TableS1, TableS1$Accession %in% unique(GW_LOH_744isolates$Sample))

#chromosome level stats were also in
ChromosomeStats <- read.csv(here("data_out", "06-Heterozygosity","GW_LOH_tables", "250813_741isolates_hetStats.csv"))

####################################
# Summary statistics
####################################

hist(TableS1_LOH$frac_HET) #744
boxplot(TableS1_LOH$frac_HET)

mean(TableS1_LOH$frac_HET) #0.006493715
sd(TableS1_LOH$frac_HET) #0.006493715
range(TableS1_LOH$frac_HET) #0.00287436 0.01390181

# by chromosome
chrStatsSum <- data.frame("chr" = paste0("chr", c(1:7, "R")), meanHet = apply(ChromosomeStats[, c("HET_ratio_chr1", "HET_ratio_chr2", "HET_ratio_chr3", "HET_ratio_chr4", "HET_ratio_chr5", "HET_ratio_chr6", "HET_ratio_chr7", "HET_ratio_chrR")], 2, function(x) mean(x)), SDhet = apply(ChromosomeStats[, c("HET_ratio_chr1", "HET_ratio_chr2", "HET_ratio_chr3", "HET_ratio_chr4", "HET_ratio_chr5", "HET_ratio_chr6", "HET_ratio_chr7", "HET_ratio_chrR")], 2, function(x) sd(x)))

write.csv(chrStatsSum, here("Calb_manu_figures_tables", "TableS2-chrStatsSummary_het.csv"), row.names=FALSE)

# by clade
LOH_clade <- split(TableS1_LOH, TableS1_LOH$Proposed_final_clade)
cladeStatsSum <- data.frame("clade" = paste0("clade", names(LOH_clade)), meanHet = unlist(lapply(LOH_clade, function(x) mean(x$frac_HET))), SDhet = unlist(lapply(LOH_clade, function(x) sd(x$frac_HET))))
range(cladeStatsSum$meanHet)

LOH_isolation <- split(TableS1_LOH, TableS1_LOH$Isolation_source)
cladeStatsSum <- data.frame("clade" = names(LOH_isolation), meanHet = unlist(lapply(LOH_isolation, function(x) mean(x$frac_HET))), SDhet = unlist(lapply(LOH_isolation, function(x) sd(x$frac_HET))))
#range(isolationStatsSum$meanHet)

order(cladeStatsSum$meanHet)


# this is from an old paper
#To examine the factors that influence genome size change, we ran two-way analyses of variance (ANOVAs) with strain, evo- lutionary environment, and their interaction as the predictor vari- ables and change in genome size as the response. When the inter- action was not significant, we reran the model without this term. To determine the relationship among environments, we ran the HSD.test function (i.e., Tukey’s test with multiple comparisons) from the agricolae package (de Mendiburu 2015).

LOHANOVA3 <- aov(frac_HET~Isolation_source*Continent*Proposed_final_clade, data = TableS1_LOH)
summary(LOHANOVA3)
# Df    Sum Sq   Mean Sq F value  Pr(>F)    
# Isolation_source                                  9 0.0001207 1.341e-05  16.508 < 2e-16 ***
#   Continent                                         4 0.0000708 1.771e-05  21.788 < 2e-16 ***
#   Proposed_final_clade                             20 0.0001855 9.275e-06  11.413 < 2e-16 ***
#   Isolation_source:Continent                       14 0.0000116 8.280e-07   1.019 0.43258    
# Isolation_source:Proposed_final_clade            85 0.0001043 1.227e-06   1.510 0.00388 ** 
#   Continent:Proposed_final_clade                   25 0.0000156 6.230e-07   0.766 0.78618    
# Isolation_source:Continent:Proposed_final_clade  19 0.0000112 5.880e-07   0.724 0.79547    
# Residuals                                       549 0.0004462 8.130e-07                    

LOHANOVA2 <- aov(frac_HET~Isolation_source*Proposed_final_clade, data = TableS1_LOH)
summary(LOHANOVA2)
# Df    Sum Sq   Mean Sq F value   Pr(>F)    
# Isolation_source                        9 0.0001278 1.420e-05  17.422  < 2e-16 ***
#   Proposed_final_clade                   20 0.0002358 1.179e-05  14.467  < 2e-16 ***
#   Isolation_source:Proposed_final_clade  85 0.0001153 1.357e-06   1.665 0.000394 ***
#   Residuals                             615 0.0005011 8.150e-07                     

out2 <- HSD.test(LOHANOVA2, "Isolation_source")
plot(out2)

##############################
# stopping here
##############################

centromere_positions <- data.frame(
  CHROM = c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chrR"),  # Chromosomes
  start = c(1563038,1927255,823334,992579,468716,980040,425812,1743190 ), # Start of centromere region
  end = c(1565967, 1930214, 826482,996216 ,471745,983792,428712,1747664)    # End of centromere region
)

centromere_positions$binStart <- round(centromere_positions$start/5000, 0)
centromere_positions$binEnd <- round(centromere_positions$end/5000, 0)

GW_LOH_744isolates_overall <- aggregate(GW_LOH_744isolates[c("Heterozygous_SNPs")], GW_LOH_744isolates[c("Bin", "CHROM")], mean, na.rm=TRUE)
hist(GW_LOH_744isolates_overall$Heterozygous_SNPs)
abline(v = mean(GW_LOH_744isolates_overall$Heterozygous_SNPs), col="red", lty=2)
abline(v = median(GW_LOH_744isolates_overall$Heterozygous_SNPs), col="blue", lty=2)
quantile(GW_LOH_744isolates_overall$Heterozygous_SNPs, probs = c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99))

#       1%       2.5%         5%        95%      97.5%        99% 
#   0.2241532  2.4431788  6.8010081 46.3916667 55.7649866 68.9421102 

#look at low and high SNP bins

#5 % high and low
mean_het_low <- subset(GW_LOH_744isolates_overall, Heterozygous_SNPs < 5)
nrow(mean_het_low)  #119 sites
table(mean_het_low$CHROM)
# chr1 chr2 chr3 chr4 chr5 chrR 
# 3    1   12    1    2   78 

mean_het_high <- subset(GW_LOH_744isolates_overall, Heterozygous_SNPs > 100)
nrow(mean_het_high) #10
table(mean_het_high$CHROM)
# chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR 
# 25   18   13   18   21   23    9   19 

write.csv(mean_het_low, here("data_out", "06-Heterozygosity", "GW_LOH_tables", "260114mean_het_low.csv"), row.names=FALSE)

write.csv(mean_het_high, here("data_out", "06-Heterozygosity", "GW_LOH_tables",  "260114mean_het_high.csv"), row.names=FALSE)

