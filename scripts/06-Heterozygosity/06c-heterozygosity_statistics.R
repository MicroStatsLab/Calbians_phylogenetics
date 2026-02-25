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

TableS1 <- read_csv("Calb_manu_tables/TableS1.csv")
GW_LOH_744isolates <- read.csv("data_out/06-Heterozygosity/250813_744isolates_GW_LOH_table.csv")

# select only the isolates that pass the HET filter
TableS1_LOH <- subset(TableS1, TableS1$Accession %in% unique(GW_LOH_744isolates$Sample))

#chromosome level stats were also in
ChromosomeStats <- read.csv(here("data_out", "06-Heterozygosity","250813_741isolates_hetStats.csv"))

####################################
# Summary statistics
####################################

hist(TableS1_LOH$frac_HET) 
boxplot(TableS1_LOH$frac_HET)

mean(TableS1_LOH$frac_HET) #0.006493715
sd(TableS1_LOH$frac_HET) #0.006493715
range(TableS1_LOH$frac_HET) #0.00287436 0.01390181

# by chromosome
chrStatsSum <- data.frame("chr" = paste0("chr", c(1:7, "R")), meanHet = apply(ChromosomeStats[, c("HET_ratio_chr1", "HET_ratio_chr2", "HET_ratio_chr3", "HET_ratio_chr4", "HET_ratio_chr5", "HET_ratio_chr6", "HET_ratio_chr7", "HET_ratio_chrR")], 2, function(x) mean(x)), SDhet = apply(ChromosomeStats[, c("HET_ratio_chr1", "HET_ratio_chr2", "HET_ratio_chr3", "HET_ratio_chr4", "HET_ratio_chr5", "HET_ratio_chr6", "HET_ratio_chr7", "HET_ratio_chrR")], 2, function(x) sd(x)))

write.csv(chrStatsSum, here("Calb_manu_tables", "TableS2-chrStatsSummary_het.csv"), row.names=FALSE)

# by clade
LOH_clade <- split(TableS1_LOH, TableS1_LOH$Proposed_final_clade)
cladeStatsSum <- data.frame("clade" = paste0("clade", names(LOH_clade)), meanHet = unlist(lapply(LOH_clade, function(x) mean(x$frac_HET))), SDhet = unlist(lapply(LOH_clade, function(x) sd(x$frac_HET))))
range(cladeStatsSum$meanHet)

LOH_isolation <- split(TableS1_LOH, TableS1_LOH$Isolation_source)
cladeStatsSum <- data.frame("clade" = names(LOH_isolation), meanHet = unlist(lapply(LOH_isolation, function(x) mean(x$frac_HET))), SDhet = unlist(lapply(LOH_isolation, function(x) sd(x$frac_HET))))
range(isolationStatsSum$meanHet)

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

GW_LOH_744isolates_mean_bin <- aggregate(GW_LOH_744isolates[c("Heterozygous_SNPs")], GW_LOH_744isolates[c("Bin", "CHROM")], mean, na.rm=TRUE)
hist(log10(GW_LOH_744isolates_mean_bin$Heterozygous_SNPs))
#look at low and high SNP bins

mean_het_low <- subset(GW_LOH_744isolates_mean_bin, Heterozygous_SNPs < 5) # 90 positions
mean_het_high <- subset(GW_LOH_744isolates_mean_bin, Heterozygous_SNPs > 100)  #27 positions
write.csv(mean_het_low, here("data_out", "06-Heterozygosity", "250815mean_het_low.csv"), row.names=FALSE)
write.csv(mean_het_high, here("data_out", "06-Heterozygosity", "250815mean_het_high.csv"), row.names=FALSE)


table(mean_het_low$CHROM)
# chr1 chr3 chr4 chr5 chrR 
# 1    9    1    1   78 

table(mean_het_high$CHROM)
# chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR 
# 120   88   69  147   65  184  113  119 

mean_het_low_all <- subset(mean_het_low, Clade == 0)
table(mean_het_low_all$CHROM)
#chrR 
#60

mean_het_high_all <- subset(mean_het_high, Clade == 0)
table(mean_het_high_all$CHROM)

# chr1 chr2 chr3 chr4 chr6 chr7 chrR 
# 4    3    1    4    7    5    3 


#did this manually for each chr x bin
which(GW_LOH_744isolates_overall$CHROM == "chr2" & GW_LOH_744isolates_overall$Bin == 386)

cenPos <- c(313, 1023, 1024, 1250, 1644, 1860, 2201, 2202, 2297, 2298, 2751, 2752)
TeloPos <- c(1, 2, 3, 638, 1085, 1445, 1766, 2005, 2212, 2402)


alldat <- data.frame(globalPos = seq_along(GW_LOH_744isolates_overall$Bin), (GW_LOH_744isolates_overall$mean_Het_SNPs))
allDatCen <- subset(alldat, globalPos %in% c(cenPos-1, cenPos, cenPos+1))
names(allDatCen)[2] <- "mean_Het_SNPs"
allDatTelo <- subset(alldat, globalPos %in% c(TeloPos-3, TeloPos-2, TeloPos-1, TeloPos, TeloPos+1, TeloPos+2))
names(allDatTelo)[2] <- "mean_Het_SNPs"

# the low cutoff is log10(5) = less than 1 SNPs per 1000 bp (5 snps in a 5kb bin)

par(mfrow=c(2, 1), mar=c(0, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(seq_along(GW_LOH_744isolates_overall$Bin), log10(GW_LOH_744isolates_overall$mean_Het_SNPs), xaxt="n", pch = ifelse(log10(GW_LOH_744isolates_overall$mean_Het_SNPs) < log10(5), 19, 21), col = ifelse(log10(GW_LOH_744isolates_overall$mean_Het_SNPs) < log10(5), "orange", "black"), yaxt="n", xlab="", ylab= "", cex=0.5, ylim=c(log10(0.01), log10(200)), xlim=c(100, 2770))
# none of the low points are on centromeres 
#points(allDatCen$globalPos, log10(allDatCen$mean_Het_SNPs), pch = ifelse(log10(allDatCen$mean_Het_SNPs) < 0.6, 19, NA), col = ifelse(log10(allDatCen$mean_Het_SNPs) < 0.6, "white", NA))
points(allDatTelo$globalPos, log10(allDatTelo$mean_Het_SNPs), pch = ifelse(log10(allDatTelo$mean_Het_SNPs) < 0.6, 19, NA), col = ifelse(log10(allDatTelo$mean_Het_SNPs) < 0.6, "white", NA))
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

plot(seq_along(GW_LOH_744isolates_overall$Bin), GW_LOH_744isolates_overall$mean_Het_SNPs, xaxt="n", pch = ifelse(GW_LOH_744isolates_overall$mean_Het_SNPs > 100, 19, 21), cex= 0.5, col = ifelse(GW_LOH_744isolates_overall$mean_Het_SNPs > 100, "purple", "black"), yaxt="n", xlab="", ylab= "", ylim=c(0, 200), xaxs="i", xlim=c(-10, 2870))
points(allDatCen$globalPos, allDatCen$mean_Het_SNPs, pch = ifelse(allDatCen$mean_Het_SNPs > 100, 19, NA), col = ifelse(allDatCen$mean_Het_SNPs > 100, "white", NA))
points(allDatTelo$globalPos, allDatTelo$mean_Het_SNPs, pch = ifelse(allDatTelo$mean_Het_SNPs > 100, 19, NA), col = ifelse(allDatTelo$mean_Het_SNPs > 100, "white", NA))
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


##### 
# Pull out high and low positions
#####

combined_het_low <- subset(GW_LOH_744isolates_overall, log10(mean_Het_SNPs) < log10(5)) # 90 positions 

# CHROM   Bin mean_Het_SNPs Clade
# <chr> <int>         <dbl> <chr>
#   1 chr1     13         3.96  0    
# 2 chr3    276         1.70  0    
# 3 chr3    278         3.89  0    
# 4 chr3    281         1.68  0    
# 5 chr3    285         3.48  0    
# 6 chr3    288         3.26  0    
# 7 chr4    136         2.02  0    
# 8 chr5    239         3.78  0   


combined_het_high <- subset(GW_LOH_744isolates_overall, mean_Het_SNPs > 100)  #27 positions

# CHROM   start     end binStart binEnd
# 1  chr1 1563038 1565967      313    313
# 2  chr2 1927255 1930214      385    386
# 3  chr3  823334  826482      165    165
# 4  chr4  992579  996216      199    199
# 5  chr5  468716  471745       94     94
# 6  chr6  980040  983792      196    197
# 7  chr7  425812  428712       85     86
# 8  chrR 1743190 1747664      349    350

write.csv(combined_het_low, "250409combined_het_low_bins.csv", row.names=FALSE)
write.csv(combined_het_high, "250409combined_het_high_bins.csv", row.names=FALSE)


#print(mean_het_high_all, n = 40)
binRanges <- GW_LOH_744isolates_overall %>% 
  group_by(CHROM) %>% 
  summarize(max = max(Bin))
# chr1 (638), 1, 337-339     
# chr2 (447), 429-432
# chr3 (360), 359
# chr4 (321), 167, 199, 240, 242
# chr5 (239), 39
# chr6 (207), 39-41, 43, 127, 160, 175
# chr7 (190), 47-48
# chrR (458), 271-272, 350

mean_het_low_all <- subset(mean_het_low, Clade == 0)

#print(mean_het_low_all, n = 60)
# chrR (458), 381-257


# PLAYING START

LOH_clade <- split(GW_LOH_744isolates, GW_LOH_744isolates$Clade)

LOH_clade1_chr6 <- subset(GW_LOH_744isolates, Clade == 1 & CHROM == "chr6")
LOH_clade1_chr6<-LOH_clade1_chr6[order(LOH_clade1_chr6$Clade_Factor),]
LOH_clade1_sample <- split(LOH_clade1_chr6, LOH_clade1_chr6$Clade_Factor)


par(mfrow=c(10, 10), mar=c(1, 1, 1, 1))
for (i in 1:100){
  print(i)
  plot(seq_along(LOH_clade1_sample[[i]]$Bin), LOH_clade1_sample[[i]]$Heterozygous_SNPs, ylim=c(0, 150), xlab="n", ylab="", cex=0.5)
  text(30, 140, i, cex=0.8)
}
par(mfrow=c(10, 10), mar=c(1, 1, 1, 1))
for (i in 101:length(LOH_clade1_sample)){
  print(i)
  plot(seq_along(LOH_clade1_sample[[i]]$Bin), LOH_clade1_sample[[i]]$Heterozygous_SNPs, ylim=c(0, 150), xlab="n", ylab="", cex=0.5)
  text(30, 140, i, cex=0.8)
}

#block 153-181
LOH_clade1_chr6_block <- subset(LOH_clade1_chr6, Clade_Factor %in% names(LOH_clade1_sample)[153:181] & Bin == 1)
LOH_clade1_chr6_block <- subset(LOH_clade1_chr6_block, Bin == 1)
table(LOH_clade1_chr6_block$Refined_isolation_source)

# Abdominal          Circulatory        Environmental 
# 1                    1                    1 
# Gastrointestinal                 Oral Pleural and Thoracic 
# 6                    6                    1 
# Skin           Urogenital               Vagina 
# 2                    4                    3 