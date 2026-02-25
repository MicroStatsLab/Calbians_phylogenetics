###############################################################################
# required libraries
###############################################################################

library(here)
library(data.table)

###############################################################################
# functions to calculate depth of coverage in bins and plot it
###############################################################################

# calculate depth of coverage in each bin (default = 5kb)
calculate_bins <- function(ddn, line, len = 5000) {
  line.ddn <- lapply(ddn, "[", c(1, 2, line))
  bin_cov <- lapply(line.ddn, function(chrom_data) {
    chrom_data$bin <- cut(chrom_data$locus, 
                          breaks = seq(0, max(chrom_data$locus), by = len), 
                          labels = FALSE)
    tapply(chrom_data[[3]], chrom_data$bin, mean, na.rm = TRUE)
  })
  names(bin_cov) <- names(ddn)
  return(bin_cov)
}

# basic graph for depth of coverage across the genome for a single line
SWline <- function(bin, title){
  #nm <-deparse(substitute(bin))
  # calculate the number of positions in the first chromosome
  x0 <- length(bin[[1]])
  # calculate the median coverage on the first chromosme 
  stand <- median(unlist(bin[[1]]))
  # plot coverage on the first chromosome, divided by the standard, multiplied by 2
  plot(1:x0,  2*unlist(bin[[1]])/stand,  type="l", col="red", 
       xlim=c(0, length(unlist(bin))), xaxt="n", yaxt="n", xlab="", ylab="", 
       ylim=c(0, 4), main = "")
  # repeat the above procedure for the remaining 7 chromosomes
  for (i in 2:8){
    x1 <- length(bin[[i]])+x0-1
    x <- x0:x1
    stand <- median(unlist(bin[[1]]))
    points(c(x0:x1), 2*unlist(bin[[i]]/stand), type="l", col=col[i])
    x0 <- x1
  }
  # titles
  mtext(title, side=3, adj=0.01, cex=0.7)  
  # plot horizontal lines at 1-4
  abline(h=1, lty=2)
  abline(h=2, lty=2)
  abline(h=3, lty=2)
  abline(h=4, lty=2)
}

# revised graph to individualize components by line (could have used this 
# for all but didn't write it until the 'closer look' section)
# enables user to change the base ploidy from two, the min and max ploidy that 
# are plotted and whether to add a y axis to the graph

SWline_closer <- function(bin, title, standFormula = "median", base = 2, minP = 0, maxP = 4, plotYaxis = FALSE){
  #nm <-deparse(substitute(bin))
  x0 <- length(bin[[1]])
  if(standFormula == "median") stand <- median(unlist(lapply(bin, function(x) median(x))))
  if(standFormula == "quant25") stand <- quantile(unlist(lapply(bin, function(x) median(x))), 0.25)
  if(standFormula == "quant75") stand <- quantile(unlist(lapply(bin, function(x) median(x))), 0.75)
  if(standFormula == "chr2") stand <- median(unlist(bin[[2]]))
  if(standFormula == "chr3") stand <- median(unlist(bin[[3]]))
  plot(1:x0,  base*unlist(bin[[1]])/stand,  type="l", col="red", 
       xlim=c(0, length(unlist(bin))), xaxt="n", yaxt="n", xlab="", ylab="", 
       ylim=c(minP, maxP), main = "")
  for (i in 2:8){
    x1 <- length(bin[[i]])+x0-1
    x <- x0:x1
    stand <- median(unlist(lapply(bin, function(x) median(x))))
    points(c(x0:x1), base*unlist(bin[[i]]/stand), type="l", col=col[i])
    x0 <- x1
  }
  mtext(title, side=3, adj=0.01, cex=0.7)  
  for(k in minP:maxP) abline(h=k, lty=3)
  if (plotYaxis) axis(2, las=2)
}

# this is required for the SWline functions to alternate red and blue chromosome colourings
col <- c(rep(c("red", "blue"), 7))


###############################################################################
# functions for processing data in batches in a loop
###############################################################################
# to load the samtools data
loaddata <- function(x) fread(here("data_in", "05-Aneuploidy", "samtools_depth", x), 
                              sep="\t", header=T)

# function to split the depth data in a loop
splitdata <- function(x){
  df <- as.data.frame(x)
  split(df, df["Chr"])
}

# to get the file names correctly
firstnum <- seq(1, 1111, by = 30)
secondnum <- seq(30, 1140, by = 30)
ranges <- paste(firstnum, secondnum, sep="_")

lookclose <- c("ERR2708450", "ERR2708456", "ERR4669765", "ERR4669777", "HSC046", "HSC050", "HSC083", "HSC122", "HSC131", "HSC139", "HSC158", "HSC160", "HSC182", "HSC210", "SRR10209617", "SRR13587535", "SRR13587542", "SRR13587569", "SRR13587577", "SRR13587579", "SRR13587596", "SRR13587612", "SRR13587639", "SRR13587643", "SRR13587657", "SRR13587673", "SRR13587680", "SRR13587748", "SRR13587782", "SRR13587761", "SRR13587868", "SRR13587878", "SRR13587889", "SRR13741117", "SRR13741118", "SRR13741121", "SRR1553997", "SRR1554177", "SRR1554178", "SRR1554287", "SRR1554289", "SRR1554296", "SRR1554297", "SRR1554304", "SRR1554310", "SRR16959048", "SRR19696103", "SRR23500659", "SRR3593469", "SRR392813", "SRR393519", "SRR397731", "SRR530262", "SRR543720", "SRR543723", "SRR543724", "SRR640897", "SRR6669863", "SRR6669867", "SRR6669868", "SRR6669891", "SRR6669909", "SRR6669920", "SRR6669922", "SRR6669940", "SRR6669942", "SRR6669948", "SRR6669989", "SRR6670027", "SRR13587749", "ERR4669752", "SRR12969478", "SRR13587712", "SRR13587806", "SRR13587818", "SRR1554302", "SRR21359696", "SRR21359704", "SRR9073715", "SRR9073716")
lookcloser2 <- c("SRR9073716", "SRR13587889", "HSC158", "SRR9073715", "SRR21359704", "ERR4669752", "SRR13587818", "SRR13587712", "SRR12969478", "SRR1554302", "SRR21359696", "SRR13587806")
finalIsol <- "SRR13587712"
###############################################################################
# read in the metadata table 
###############################################################################
metadata <- read.csv(here("Calb_manu_tables", "TableS1.csv"), header = TRUE) #1180

# subset metadata table to get names of isolates that passed the required filter for this analysis
metadata_pass <- subset(metadata, PASS_Het == 1) #842 #250907 908

metadata_aneup <- subset(metadata, Accession %in% lookclose) #842
write.csv(metadata_aneup, here("Calb_manu_tables", "250626TableSx_karyotypicVars.csv")) #1180)

initcol <-  fread(here("data_in", "05-Aneuploidy","depth_1_30_coverage.txt"), sep="\t", header=T)

finalcol <-  fread(here("data_in", "05-Aneuploidy","samtools_depth", "coverage_batch2", "depth_301_330_coverage.txt"), sep="\t", header=T)

finalfinalcol <-  fread(here("data_in", "05-Aneuploidy","samtools_depth", "250806_redo_sample_depth.txt"), sep="\t", header=T)


###############################################################################
# read in look closerlook data
###############################################################################
lookclose_data_1_270 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_1_270.csv"), sep=",", header=T)  
lookclose_data_1_270 <- cbind(initcol[,1:2], lookclose_data_1_270)
lookclose_data_1_270 <- as.data.frame(lookclose_data_1_270)

lookclose_data_271_570 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_271_570.csv"), sep=",", header=T)  
lookclose_data_271_270 <- cbind(initcol[,1:2], lookclose_data_271_570)
lookclose_data_271_270 <- as.data.frame(lookclose_data_271_270)

lookclose_data_571_870 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_571_870.csv"), sep=",", header=T)  
lookclose_data_571_870 <- cbind(initcol[,1:2], lookclose_data_571_870)
lookclose_data_571_870 <- as.data.frame(lookclose_data_571_870)

lookclose_data_871_1140 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_871_1140.csv"), sep=",", header=T)  
lookclose_data_871_1140 <- cbind(initcol[,1:2], lookclose_data_871_1140)
lookclose_data_871_1140 <- as.data.frame(lookclose_data_871_1140)

###############################################################################
# calculate bins again from closer look data
# set one, isol 1-270
###############################################################################
isol_1_270 <- names(lookclose_data_1_270)[4:length(lookclose_data_1_270)]
for(i in 1:length(isol_1_270)){
  splitname <- paste0(isol_1_270[i], "_ddn")
  #print(splitname)
  dataset <- lookclose_data_1_270[, isol_1_270[i]]
  dataset <- cbind(lookclose_data_1_270[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_1_270[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}


###############################################################################
# generate graphs to look at closerlook depths
# for 'complicated' karyotypes, conduct individual assessment 
# change base ploidy as needed, play with bin size, etc.

# NB - this is the only place where individual isolate names are used, 
# so be extra careful not to make mistakes here
###############################################################################
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook1_20.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
SWline_closer(ERR2708450_bins, title = "ERR2708450") # just messy/U can't score
SWline_closer(ERR2708456_bins, title = "ERR2708456") # just messy/U can't score
SWline_closer(ERR4669765_bins, title = "ERR4669765") # +chr7
SWline_closer(ERR4669777_bins, title = "ERR4669777") # +chrR- right arm, 1362000
SWline_closer(HSC046_bins, title = "HSC046") # +chr7
SWline_closer(HSC050_bins, title = "HSC050") # +chr7
SWline_closer(HSC083_bins, title = "HSC083") # chr4 + chr7
SWline_closer(HSC122_bins, title = "HSC122") # +chr1 almost all +chr5 L
SWline_closer(HSC131_bins, title = "HSC131") # +chr6 +chr7 most
SWline_closer(HSC139_bins, title = "HSC139") # + chr6x2 most +chr7x2

SWline_closer(HSC158_bins, title = "HSC158") # +chrR right
SWline_closer(HSC160_bins, title = "HSC160") # +chr6 +chr7 -> +chr6 mixed pop? Nothing clean
SWline_closer(HSC182_bins, title = "HSC182") # +chr7
SWline_closer(HSC210_bins, title = "HSC210") # 4N -chr7 
SWline_closer(SRR10209617_bins, title = "SRR10209617") # +chr2 slightly elevated? mixed?
SWline_closer(SRR13587535_bins, title = "SRR13587535") # 4N +chr7
SWline_closer(SRR13587542_bins, title = "SRR13587542") # -chr5 part, -chr6 part
SWline_closer(SRR13587569_bins, title = "SRR13587569") # +chrR
SWline_closer(SRR13587577_bins, title = "SRR13587577") # +chr7 left part +chrR right
SWline_closer(SRR13587579_bins, title = "SRR13587579") # +chr6 +chr7 - not diploid?
dev.off()

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook21_24.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
# set 1, 1-270
SWline_closer(SRR13587596_bins, title = "SRR13587596") # +chr7
SWline_closer(SRR13587612_bins, title = "SRR13587612") # +chr7
SWline_closer(SRR13587639_bins, title = "SRR13587639") # +chr4 part +chrR part
SWline_closer(SRR13587643_bins, title = "SRR13587643") # +chr2
dev.off()


# ERR2708450_bins
SWline_closer(ERR2708450_bins, title = "ERR2708450") # try different stand

ERR2708450_bins_10 <- calculate_bins(ERR2708450_ddn, 3, len = 10000)
medianDepth <- lapply(ERR2708450_bins_10, function(x) median(x))
quant25Depth <- lapply(ERR2708450_bins_10, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "ERR2708450_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(1, 2), mar=c(3, 3, 3, 0))
plot(1:8, unlist(medianDepth)/min(unlist(medianDepth))*2, ylim=c(1, 6), pch=19, main = "median depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(medianDepth)/min(unlist(medianDepth))*3, pch=19, col=grey(0.4))
points(unlist(medianDepth)/min(unlist(medianDepth))*4, pch=19, col=grey(0.8))

plot(1:8, unlist(quant25Depth)/min(unlist(quant25Depth))*2, ylim=c(1, 6), pch=19, main = "25 quantile depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(quant25Depth)/min(unlist(quant25Depth))*3, pch=19, col=grey(0.4))
points(unlist(quant25Depth)/min(unlist(quant25Depth))*4, pch=19, col=grey(0.8))
dev.off()

# ERR2708456_bins
SWline_closer(ERR2708456_bins, title = "ERR2708456") # try different stand

ERR2708456_bins_10 <- calculate_bins(ERR2708456_ddn, 3, len = 10000)
medianDepth <- lapply(ERR2708456_bins_10, function(x) median(x))
quant25Depth <- lapply(ERR2708456_bins_10, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "ERR2708456_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(1, 2), mar=c(3, 3, 3, 0))
plot(1:8, unlist(medianDepth)/min(unlist(medianDepth))*2, ylim=c(1, 6), pch=19, main = "median depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(medianDepth)/min(unlist(medianDepth))*3, pch=19, col=grey(0.4))
points(unlist(medianDepth)/min(unlist(medianDepth))*4, pch=19, col=grey(0.8))

plot(1:8, unlist(quant25Depth)/min(unlist(quant25Depth))*2, ylim=c(1, 6), pch=19, main = "25 quantile depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(quant25Depth)/min(unlist(quant25Depth))*3, pch=19, col=grey(0.4))
points(unlist(quant25Depth)/min(unlist(quant25Depth))*4, pch=19, col=grey(0.8))
dev.off()

# ERR4669777_bins
SWline_closer(ERR4669777_bins, title = "ERR4669777") # try different stand
quant25Depth <- lapply(ERR4669777_bins, function(x) quantile(x, 0.25)) #for partial aneuploidies - use 25% quantile for the standardization

ERR4669777_bins_25 <- calculate_bins(ERR4669777_ddn, 3, len = 2500)
ERR4669777_bins_20 <- calculate_bins(ERR4669777_ddn, 3, len = 2000)

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "ERR4669777_depth.pdf"), width=4, height=6, font = "Times")
plot(ERR4669777_bins_20[[8]]/quant25Depth[[8]]*2, ylim=c(0, 4))
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(ERR4669777_bins_20[[8]]/quant25Depth[[8]]*2 > 2.75)
abline(v=688, col="red") # approx 1362000
dev.off()

# HSC122
SWline_closer(HSC122_bins, title = "HSC122") # try different stand
quant25Depth <- lapply(HSC122_bins, function(x) quantile(x, 0.25)) #for partial aneuploidies - use 25% quantile for the standardization
quant75Depth <- lapply(HSC122_bins, function(x) quantile(x, 0.75)) #for partial aneuploidies - use 25% quantile for the standardization

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC122_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(HSC122_bins[[1]]/quant75Depth[[1]]*2, ylim=c(0, 4), main = "Chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(HSC122_bins[[1]]/quant75Depth[[1]]*2 < 1.5)
abline(v=92, col="red") # approx 460000

plot(HSC122_bins[[5]]/quant75Depth[[2]]*2, ylim=c(0, 5), main = "Chr5")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
which(HSC122_bins[[5]]/quant75Depth[[2]]*2 >3)
abline(v=94, col="red") # approx 470000
dev.off()

# HSC131
SWline_closer(HSC131_bins, title = "HSC131") # try different stand
quant25Depth <- lapply(HSC131_bins, function(x) quantile(x, 0.25))
quant50Depth <- lapply(HSC131_bins, function(x) quantile(x, 0.5)) 
quant75Depth <- lapply(HSC131_bins, function(x) quantile(x, 0.75))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC131_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(HSC131_bins[[6]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC131_bins[[6]]/quant50Depth[[1]]*2 < 2.5)
abline(v=173, col="red") # approx 865000

plot(HSC131_bins[[7]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC131_bins[[7]]/quant50Depth[[1]]*2 > 2.5)
abline(v=46, col="red") # approx 230000
dev.off()

# HSC139
SWline_closer(HSC139_bins, title = "HSC139") # try different stand
quant25Depth <- lapply(HSC139_bins, function(x) quantile(x, 0.25))
quant50Depth <- lapply(HSC139_bins, function(x) quantile(x, 0.5)) 
quant75Depth <- lapply(HSC139_bins, function(x) quantile(x, 0.75))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC139_depth.pdf"), width=4, height=8, font = "Times")
par(mfrow=c(3, 1))
plot(HSC139_bins[[5]]/quant50Depth[[1]]*2, ylim=c(0, 6), main = "Chr5")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
which(HSC139_bins[[5]]/quant50Depth[[5]]*2 >5)
abline(v=94, col="red") # approx 470000

plot(HSC139_bins[[6]]/quant50Depth[[1]]*2, ylim=c(0, 4.5), main = "Chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC139_bins[[6]]/quant50Depth[[1]]*2 > 2.5)
abline(v=41, col="red") # 205000 
abline(v=175, col="red") # 875000

plot(HSC139_bins[[7]]/quant50Depth[[1]]*2, ylim=c(0, 5), main = "Chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
#which(HSC139_bins[[7]]/quant50Depth[[1]]*2 > 2.5)
#abline(v=41, col="red") # 205000 
dev.off()

# HSC158
SWline_closer(HSC158_bins, title = "HSC158") # try different stand
quant50Depth <- lapply(HSC158_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC158_depth.pdf"), width=4, height=6, font = "Times")
plot(HSC158_bins[[8]]/quant50Depth[[2]]*2, ylim=c(0, 4))
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(HSC158_bins[[8]]/quant50Depth[[2]]*2 > 2.75)
abline(v=226, col="red") # approx 1130000
dev.off()


# HSC160
SWline_closer(HSC160_bins, title = "HSC160") # try different stand
quant50Depth <- lapply(HSC160_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "HSC160_varPloidy.pdf"), width=4, height=8, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 3, 1, 1))
SWline_closer(HSC160_bins, "HSC160, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

# HSC210
SWline_closer(HSC210_bins, title = "HSC210") # try different stand
quant50Depth <- lapply(HSC210_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC210_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 3, 1, 1))
SWline_closer(HSC210_bins, "HSC210, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(HSC210_bins, "HSC210, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(HSC210_bins, "HSC210, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()


# SRR10209617
SWline_closer(SRR10209617_bins, title = "SRR10209617") # try different stand
quant50Depth <- lapply(SRR10209617_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR10209617_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 3, 1, 1))
SWline_closer(SRR10209617_bins, "SRR10209617, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR10209617_bins, "SRR10209617, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR10209617_bins, "SRR10209617, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

# SRR13587535
SWline_closer(SRR13587535_bins, title = "SRR13587535") # try different stand
quant50Depth <- lapply(SRR13587535_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587535-depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 3, 1, 1))
SWline_closer(SRR13587535_bins, "SRR13587535, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587535_bins, "SRR13587535, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587535_bins, "SRR13587535, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

# SRR13587542
SWline_closer(SRR13587542_bins, title = "SRR13587542") # try different stand
quant50Depth <- lapply(SRR13587542_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587542_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR13587542_bins[[5]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr5")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587542_bins[[5]]/quant50Depth[[1]]*2 < 1.5)
abline(v=147, col="red") # approx 735000

plot(SRR13587542_bins[[6]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587542_bins[[6]]/quant50Depth[[1]]*2 < 1.5)
abline(v=176, col="red") # approx 880000
dev.off()

# SRR13587577
SWline_closer(SRR13587577_bins, title = "SRR13587577") # try different stand
quant50Depth <- lapply(SRR13587577_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR135875772_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR13587577_bins[[7]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587577_bins[[7]]/quant50Depth[[1]]*2 > 2.5)
abline(v=49, col="red") # approx 245000

plot(SRR13587577_bins[[8]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "ChrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587577_bins[[8]]/quant50Depth[[1]]*2 > 2.5)
abline(v=271, col="red") # approx 1355000
dev.off()

# SRR13587579_bins
SWline_closer(SRR13587579_bins, title = "SRR13587579") # try different stand
quant50Depth <- lapply(SRR13587579_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587579-depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 3, 1, 1))
SWline_closer(SRR13587579_bins, "SRR13587579, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587579_bins, "SRR13587579, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587579_bins, "SRR13587579, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

# SRR13587639
SWline_closer(SRR13587639_bins, title = "SRR13587639") # try different stand
quant50Depth <- lapply(SRR13587639_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587639_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR13587639_bins[[4]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "Chr4")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587639_bins[[4]]/quant50Depth[[1]]*2 > 2.5)
abline(v=237, col="red") # approx 1185000

plot(SRR13587639_bins[[8]]/quant50Depth[[1]]*2, ylim=c(0, 4), main = "ChrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587639_bins[[8]]/quant50Depth[[1]]*2 > 2.5)
abline(v=271, col="red") # approx 1355000
dev.off()

###############################################################################
# calculate bins again from closer look data
# set two, isol 271-570
###############################################################################

isol_271_570 <- names(lookclose_data_271_270)[4:length(lookclose_data_271_270)]
for(i in 1:length(isol_271_570)){
  splitname <- paste0(isol_271_570[i], "_ddn")
  #print(splitname)
  dataset <- lookclose_data_271_270[, isol_271_570[i]]
  dataset <- cbind(lookclose_data_271_270[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_271_570[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}

###############################################################################
# generate graphs to look at closerlook depths
# for 'complicated' karyotypes, conduct individual assessment 
# change base ploidy as needed, play with bin size, etc.

# NB - this is the only place where individual isolate names are used, 
# so be extra careful not to make mistakes here
###############################################################################
isol_271_570

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook25_47.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(11, 2), mar=c(1, 1, 1, 1))
SWline_closer(SRR13587657_bins, title = "SRR13587657") # 3N +chr1 +chr2 +chr7 (3N?)
SWline_closer(SRR13587673_bins, title = "SRR13587673") # +chr7
SWline_closer(SRR13587680_bins, title = "SRR13587680") # +chr2 +chr3 +chr6 +chrR (???)
SWline_closer(SRR13587748_bins, title = "SRR13587748") # +chr4
SWline_closer(SRR13587761_bins, title = "SRR13587761") # +chr7
SWline_closer(SRR13587782_bins, title = "SRR13587782") # +chr7
SWline_closer(SRR13587868_bins, title = "SRR13587868") # +chr1x3 +chr5x2 +chr6x2 +chr7x3 (3N)
SWline_closer(SRR13587878_bins, title = "SRR13587878") # +chr5 +chr7
SWline_closer(SRR13587889_bins, title = "SRR13587889") # +chr4 R  
SWline_closer(SRR13741117_bins, title = "SRR13741117") # +complicated + U impossible

SWline_closer(SRR13741118_bins, title = "SRR13741118") # +chr4 +chrR
SWline_closer(SRR13741121_bins, title = "SRR13741121") # U
SWline_closer(SRR1553997_bins, title = "SRR1553997") # 3N -chr4 -chr6 -chr7
SWline_closer(SRR1554177_bins, title = "SRR1554177") # +chr5 U
SWline_closer(SRR1554178_bins, title = "SRR1554178") # +chr5
SWline_closer(SRR1554287_bins, title = "SRR1554287") # +chr4 +chr5
SWline_closer(SRR1554289_bins, title = "SRR1554289") # +chr4
SWline_closer(SRR1554296_bins, title = "SRR1554296") # +chr5
SWline_closer(SRR1554297_bins, title = "SRR1554297") # +chr5
SWline_closer(SRR1554304_bins, title = "SRR1554304") # +chr5R  
SWline_closer(SRR1554310_bins, title = "SRR1554310") # +chr2R  +chr4 +chr5L
SWline_closer(SRR16959048_bins, title = "SRR16959048") #U
dev.off()

# SRR13587657
SWline_closer(SRR13587657_bins, title = "SRR13587657") # try different stand
SRR13587657_bins_10 <- calculate_bins(SRR13587657_ddn, 3, len = 10000)
medianDepth <- lapply(SRR13587657_bins_10, function(x) median(x))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR13587657_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(3, 3, 3, 0))
SWline_closer(SRR13587657_bins_10, "SRR13587657, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587657_bins_10, "SRR13587657, base = 3", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587657_bins_10, "SRR13587657, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

# SRR13587680
SWline_closer(SRR13587680_bins, title = "SRR13587680") # try different stand
SRR13587680_bins_10 <- calculate_bins(SRR13587680_ddn, 3, len = 10000)

medianDepth <- lapply(SRR13587680_bins, function(x) median(x))
#quant25Depth <- lapply(ERR2708456_bins_10, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR13587680_depth.pdf"), width=8, height=6, font = "Times")
par(mfrow=c(3, 3), mar=c(3, 3, 3, 0))
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 2, quant25", standFormula = "quant25", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 3, quant25", standFormula = "quant25", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 4, quant25", standFormula = "quant25", base = 4, 2, 6, plotYaxis = TRUE)

SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 2, quant50", standFormula = "median",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 3, quant50", standFormula = "median",  base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 4, quant50", standFormula = "median",  base = 4, 2, 6, plotYaxis = TRUE)

SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 2, quant75", standFormula = "quant75", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 3, quant75", standFormula = "quant75", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587680_bins_10, "SRR13587680, base = 4, quant75", standFormula = "quant75", base = 4, 2, 6, plotYaxis = TRUE)
dev.off()

# SRR13587868
SWline_closer(SRR13587868_bins, title = "SRR13587868") # try different stand
SRR13587868_bins_10 <- calculate_bins(SRR13587868_ddn, 3, len = 10000)

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR13587868_depth.pdf"), width=8, height=6, font = "Times")
par(mfrow=c(3, 3), mar=c(3, 3, 3, 0))
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 2, quant25", standFormula = "quant25",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 3, quant25", standFormula = "quant25",  base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 4, quant25", standFormula = "quant25",  base = 4, 2, 7, plotYaxis = TRUE)

SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 2, quant50", standFormula = "median",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 3, quant50", standFormula = "median",  base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 4, quant50", standFormula = "median",  base = 4, 2, 7, plotYaxis = TRUE)

SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 2, chr2", standFormula = "chr2",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 3, chr2", standFormula = "chr2",  base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13587868_bins_10, "SRR13587868, base = 4, chr2", standFormula = "chr2",  base = 4, 2, 7, plotYaxis = TRUE)
dev.off()

# SRR13741117
SWline_closer(SRR13741117_bins, title = "SRR13741117") # try different stand
SRR13741117_bins_10 <- calculate_bins(SRR13741117_ddn, 3, len = 10000)

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR13741117_depth.pdf"), width=8, height=6, font = "Times")
par(mfrow=c(3, 3), mar=c(3, 3, 3, 0))
SWline_closer(SRR13741117_bins_10, "SRR13741117, base = 2, quant25", standFormula = "quant25",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR13741117_bins_10, "SRR13741117, base = 3, quant25", standFormula = "quant25",  base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR13741117_bins_10, "SRR13741117, base = 4, quant25", standFormula = "quant25",  base = 4, 2, 7, plotYaxis = TRUE)
dev.off()

# SRR1553997
SWline_closer(SRR1553997_bins, title = "SRR1553997") # -chr5 -chr6 -chr7

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR1553997_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(3, 3, 3, 0))
SWline_closer(SRR1553997_bins, "SRR1553997, base = 2",standFormula = "median",  base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(SRR1553997_bins, "SRR1553997, base = 3", standFormula = "median", base = 3, 1, 6, plotYaxis = TRUE)
SWline_closer(SRR1553997_bins, "SRR1553997, base = 4", standFormula = "median", base = 4, 2, 8, plotYaxis = TRUE)
dev.off()

#SRR1554177
SWline_closer(SRR1554177_bins, title = "SRR1554177") # +chr5 U

medianDepth <- lapply(SRR1554177_bins, function(x) median(x))
quant25Depth <- lapply(SRR1554177_bins, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR1554177_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(1, 2), mar=c(3, 3, 3, 0))
plot(1:8, unlist(medianDepth)/min(unlist(medianDepth))*2, ylim=c(1, 6), pch=19, main = "median depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(medianDepth)/min(unlist(medianDepth))*3, pch=19, col=grey(0.4))
points(unlist(medianDepth)/min(unlist(medianDepth))*4, pch=19, col=grey(0.8))

plot(1:8, unlist(quant25Depth)/min(unlist(quant25Depth))*2, ylim=c(1, 6), pch=19, main = "25 quantile depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(quant25Depth)/min(unlist(quant25Depth))*3, pch=19, col=grey(0.4))
points(unlist(quant25Depth)/min(unlist(quant25Depth))*4, pch=19, col=grey(0.8))
dev.off()

# SRR1554304
SWline_closer(SRR1554304_bins, title = "SRR1554304") 
quant25Depth <- lapply(SRR1554304_bins, function(x) quantile(x, 0.25)) #for partial aneuploidies - use 25% quantile for the standardization

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR1554304_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR1554304_bins[[5]]/quant25Depth[[2]]*2, ylim=c(0, 5))
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(SRR1554304_bins[[5]]/quant25Depth[[2]]*2 > 2.75)
abline(v=95, col="red") # approx 475000
dev.off()

# SRR1554310
SWline_closer(SRR1554310_bins, title = "SRR1554310") #
medDepth <- lapply(SRR1554310_bins, function(x) quantile(x, 0.5)) #for partial aneuploidies - use 25% quantile for the standardization

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR1554310_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR1554310_bins[[2]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chr2")
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR1554310_bins[[2]]/medDepth[[3]]*2 > 3.75)
abline(v=375, col="red") # approx 475000

plot(SRR1554310_bins[[5]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chr5")
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR1554310_bins[[5]]/medDepth[[3]]*2 > 3.75)
abline(v=78, col="red") # approx 390000
dev.off()

###############################################################################
# calculate bins again from closer look data
# set three, isol 571 - 870
###############################################################################

isol_571_870 <- names(lookclose_data_571_870)[3:length(lookclose_data_571_870)]
for(i in 1:length(isol_571_870)){
  splitname <- paste0(isol_571_870[i], "_ddn")
  #print(splitname)
  dataset <- lookclose_data_571_870[, isol_571_870[i]]
  dataset <- cbind(lookclose_data_571_870[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_571_870[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook48_49.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
SWline_closer(SRR19696103_bins, title = "SRR19696103") #-chr6R +chrR-R
SWline_closer(SRR23500659_bins, title = "SRR23500659") # mess
dev.off()

#SRR19696103
SWline_closer(SRR19696103_bins, title = "SRR19696103") 
medDepth <- lapply(SRR19696103_bins, function(x) quantile(x, 0.5)) #for partial aneuploidies - use 25% quantile for the standardization

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR19696103_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR19696103_bins[[6]]/medDepth[[2]]*2, ylim=c(0, 5), main = "chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(SRR19696103_bins[[6]]/medDepth[[2]]*2 < 1.5)
abline(v=176, col="red") # approx 880000

plot(SRR19696103_bins[[8]]/medDepth[[2]]*2, ylim=c(0, 5), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(SRR19696103_bins[[8]]/medDepth[[2]]*2 > 2.75)
abline(v=287, col="red") # approx 1435000
dev.off()


#SRR23500659
SWline_closer(SRR23500659_bins, title = "SRR23500659") 

medianDepth <- lapply(SRR23500659_bins, function(x) median(x))
quant25Depth <- lapply(SRR23500659_bins, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", "SRR23500659_depth.pdf"), width=6, height=6, font = "Times")
par(mfrow=c(1, 2), mar=c(3, 3, 3, 0))
plot(1:8, unlist(medianDepth)/min(unlist(medianDepth))*2, ylim=c(1, 6), pch=19, main = "median depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(medianDepth)/min(unlist(medianDepth))*3, pch=19, col=grey(0.4))
points(unlist(medianDepth)/min(unlist(medianDepth))*4, pch=19, col=grey(0.8))

plot(1:8, unlist(quant25Depth)/min(unlist(quant25Depth))*2, ylim=c(1, 6), pch=19, main = "25 quantile depth", yaxt="n")
axis(2, las=2)
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)
abline(h=6, lty=3)
points(unlist(quant25Depth)/min(unlist(quant25Depth))*3, pch=19, col=grey(0.4))
points(unlist(quant25Depth)/min(unlist(quant25Depth))*4, pch=19, col=grey(0.8))
dev.off()

###############################################################################
# calculate bins again from closer look data
# set three, isol 871 - 1140
###############################################################################

isol_871_1140 <- names(lookclose_data_871_1140)[3:length(lookclose_data_871_1140)]
for(i in 1:length(isol_871_1140)){
  splitname <- paste0(isol_871_1140[i], "_ddn")
  #print(splitname)
  dataset <- lookclose_data_871_1140[, isol_871_1140[i]]
  dataset <- cbind(lookclose_data_871_1140[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_871_1140[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}

isol_871_1140

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook50_72.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(11, 2), mar=c(1, 1, 1, 1))
SWline_closer(SRR3593469_bins, title = "SRR3593469") # U
SWline_closer(SRR392813_bins, title = "SRR392813") # +chr4
SWline_closer(SRR393519_bins, title = "SRR393519") # +chr5
SWline_closer(SRR397731_bins, title = "SRR397731") # +chr7
SWline_closer(SRR530262_bins, title = "SRR530262") # +chr5 +chr7 +complicated chrR-R
SWline_closer(SRR543720_bins, title = "SRR543720") # +chr7
SWline_closer(SRR543723_bins, title = "SRR543723") # +chr1 part
SWline_closer(SRR543724_bins, title = "SRR543724") # +chr4
SWline_closer(SRR640897_bins, title = "SRR640897") # +chr6R +chr7R
SWline_closer(SRR6669863_bins, title = "SRR6669863") # +chrRR

SWline_closer(SRR6669867_bins, title = "SRR6669867") # complicated AF
SWline_closer(SRR6669868_bins, title = "SRR6669868") # complicated AF
SWline_closer(SRR6669891_bins, title = "SRR6669891") # +chr5 +chrR
SWline_closer(SRR6669909_bins, title = "SRR6669909") # +chr7
SWline_closer(SRR6669920_bins, title = "SRR6669920") # +chr1 small +chr2
SWline_closer(SRR6669922_bins, title = "SRR6669922") # 4N +chrR 
SWline_closer(SRR6669940_bins, title = "SRR6669940") # +chr5
SWline_closer(SRR6669942_bins, title = "SRR6669942") # +chr5 not 2N
SWline_closer(SRR6669948_bins, title = "SRR6669948") # +chr7
SWline_closer(SRR6669989_bins, title = "SRR6669989") # +chr6
SWline_closer(SRR6670027_bins, title = "SRR6670027") # +chr1
SWline_closer(SRR13587749_bins, title = "SRR13587749") # +chr4
dev.off()

#SRR530262
SWline_closer(SRR530262_bins, title = "SRR530262") # +chr5 +chr7 +complicated chrR-R
medDepth <- lapply(SRR530262_bins, function(x) quantile(x, 0.5)) 

SRR530262_bins_10 <- calculate_bins(SRR530262_ddn, 3, len = 10000)
SWline_closer(SRR530262_bins_10, title = "SRR530262", base = 2, maxP = 6, ) # +chr5 +chr7 +complicated chrR-R

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR530262_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR530262_bins_10[[8]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR530262_bins_10[[8]]/medDepth[[2]]*2 > 3.5)
abline(v=197, col="red") # approx 985000
which(SRR530262_bins_10[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=149, col="red") # approx 745000
dev.off()

#SRR543723
SWline_closer(SRR543723_bins, title = "SRR543723") # +chr1 part
medDepth <- lapply(SRR543723_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR543723_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR543723_bins[[1]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR543723_bins[[1]]/medDepth[[2]]*2 > 2.5)
abline(v=87, col="red") # approx 435000
which(SRR530262_bins_10[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=170, col="red") # approx 850000
dev.off()

#SRR640897
SWline_closer(SRR640897_bins, title = "SRR640897") # +chr6R +chr7R
medDepth <- lapply(SRR640897_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR640897_depth.pdf"), width=4, height=6, font = "Times")
par(mfrow=c(2, 1))
plot(SRR640897_bins[[6]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR640897_bins[[6]]/medDepth[[2]]*2 > 2.5)
abline(v=174, col="red") # approx 435000

plot(SRR640897_bins[[7]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR640897_bins[[7]]/medDepth[[2]]*2 > 2.5)
abline(v=117, col="red") # approx 585000
dev.off()

#SRR6669863
SWline_closer(SRR6669863_bins, title = "SRR6669863") # +chrRR
medDepth <- lapply(SRR6669863_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669863_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR6669863_bins[[8]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669863_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=232, col="red") # approx 1160000
dev.off()

#SRR6669867
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669867_depth.pdf"), width=8, height=10, font = "Times")
par(mfrow=c(3, 2))
SWline_closer(SRR6669867_bins, title = "SRR6669867", base = 3, maxP = 6, plotYaxis = TRUE) # complicated AF
medDepth <- lapply(SRR6669867_bins, function(x) quantile(x, 0.5)) 

plot(SRR6669867_bins[[1]]/medDepth[[5]]*3, ylim=c(0, 5), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669867_bins[[1]]/medDepth[[5]]*3 > 2.5)
abline(v=96, col="red") # approx 480000
abline(v=605, col="red") # approx 3025000

plot(SRR6669867_bins[[4]]/medDepth[[5]]*3, ylim=c(0, 5), main = "chr4")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669867_bins[[4]]/medDepth[[5]]*3 > 2.5)
abline(v=82, col="red") # approx 410000
abline(v=170, col="red") # approx 850000
abline(v=266, col="red") # approx 1330000

plot(SRR6669867_bins[[6]]/medDepth[[5]]*3, ylim=c(0, 5), main = "chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669867_bins[[6]]/medDepth[[5]]*3 > 2.5)
abline(v=66, col="red") # approx 330000
abline(v=175, col="red") # approx 875000

plot(SRR6669867_bins[[7]]/medDepth[[5]]*3, ylim=c(0, 5), main = "chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669867_bins[[7]]/medDepth[[5]]*3 > 2.5)
abline(v=46, col="red") # approx 230000

plot(SRR6669867_bins[[8]]/medDepth[[5]]*3, ylim=c(0, 5), main = "chr8")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669867_bins[[8]]/medDepth[[5]]*3 > 2.5)
abline(v=47, col="red") # approx 235000
abline(v=97, col="red") # approx 485000
abline(v=107, col="red") # approx 535000
dev.off()

#SRR6669868
SWline_closer(SRR6669868_bins, title = "SRR6669868") # complicated AF
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669868_depth_ploidy.pdf"), width=8, height=10, font = "Times")
par(mfrow=c(3, 1))
SWline_closer(SRR6669868_bins, title = "SRR6669868", base = 2, maxP = 7, plotYaxis = TRUE, standFormula = "chr2") # complicated AF
SWline_closer(SRR6669868_bins, title = "SRR6669868", base = 3, maxP = 7, plotYaxis = TRUE, standFormula = "chr2") # complicated AF
SWline_closer(SRR6669868_bins, title = "SRR6669868", base = 4, maxP = 7, plotYaxis = TRUE, standFormula = "chr2") # complicated AF
dev.off()

medDepth <- lapply(SRR6669868_bins, function(x) quantile(x, 0.5)) 
quant25Depth <- lapply(SRR6669868_bins, function(x) quantile(x, 0.25))

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669868_depth_aneup.pdf"), width=8, height=10, font = "Times")
par(mfrow=c(3, 2))
SWline_closer(SRR6669868_bins, title = "SRR6669868", base = 4, maxP = 7, plotYaxis = TRUE, standFormula = "chr2") # complicated AF

plot(SRR6669868_bins[[1]]/medDepth[[2]]*4, ylim=c(0, 7), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
which(SRR6669868_bins[[1]]/quant25Depth[[2]]*4 > 4.5)
abline(v=239, col="red") # approx 1195000

plot(SRR6669868_bins[[4]]/quant25Depth[[2]]*4, ylim=c(0, 7), main = "chr4")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
which(SRR6669868_bins[[4]]/quant25Depth[[2]]*4 > 4.5)
abline(v=242, col="red") # approx 1210000

plot(SRR6669868_bins[[5]]/quant25Depth[[2]]*4, ylim=c(0, 8), main = "chr5")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
abline(h = 7, lty=3)
abline(h = 8, lty=3)
which(SRR6669868_bins[[5]]/quant25Depth[[2]]*4 > 7.5)
abline(v=94, col="red") # approx 470000

plot(SRR6669868_bins[[7]]/quant25Depth[[2]]*4, ylim=c(0, 8), main = "chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
which(SRR6669868_bins[[7]]/quant25Depth[[2]]*4 > 5.5)
abline(v=50, col="red") # approx 250000

plot(SRR6669868_bins[[8]]/quant25Depth[[2]]*4, ylim=c(0, 8), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
abline(h = 5, lty=3)
abline(h = 6, lty=3)
which(SRR6669868_bins[[8]]/quant25Depth[[2]]*4 > 5.5)
abline(v=15, col="red") # approx 75000
abline(v=353, col="red") # approx 1765000
dev.off()

#SRR6669920
SWline_closer(SRR6669920_bins, title = "SRR6669920") # +chr1 small +chr2
medDepth <- lapply(SRR6669920_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669920_depth.pdf"), width=8, height=10, font = "Times")
plot(SRR6669920_bins[[1]]/medDepth[[3]]*2, ylim=c(0, 4), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
which(SRR6669920_bins[[1]]/medDepth[[3]]*2 > 2.5)
abline(v=81, col="red") # approx 405000
abline(v=92, col="red") # approx 460000
abline(v=359, col="red") # approx 1795000
abline(v=379, col="red") # approx 1895000
dev.off()

# SRR6669922
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669922_depth_ploidy.pdf"), width=5, height=6, font = "Times")
par(mfrow=c(3, 1))
SWline_closer(SRR6669922_bins, title = "SRR6669922", base = 2, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR6669922_bins, title = "SRR6669922", base = 3, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR6669922_bins, title = "SRR6669922", base = 4, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
dev.off()

SWline_closer(SRR6669922_bins, title = "SRR6669922", base = 3) # +chrR not 2N

###############################################################################
# read in look closerlook2 data
###############################################################################
lookclose2_data_1_270 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose2_data_1_270.csv"), sep=",", header=T)  
lookclose2_data_271_570 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose2_data_271_570.csv"), sep=",", header=T)  
lookclose2_data_571_870 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_571_870.csv"), sep=",", header=T)  
lookclose2_data_871_1140 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose2_data_871_1140.csv"), sep=",", header=T)  

lookclose2_all_data <- cbind(initcol[,1:2], lookclose2_data_1_270, lookclose2_data_271_570, lookclose2_data_571_870, lookclose2_data_871_1140)
lookclose2_all_data <- as.data.frame(lookclose2_all_data)

isol_lookclose2 <- names(lookclose2_all_data)[3:length(lookclose2_all_data)]
for(i in 1:length(isol_lookclose2)){
  splitname <- paste0(isol_lookclose2[i], "_ddn")
  dataset <- lookclose2_all_data[, isol_lookclose2[i]]
  dataset <- cbind(lookclose2_all_data[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_lookclose2[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}


pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook2.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(6, 2), mar=c(1, 1, 1, 1))
SWline_closer(ERR4669752_bins, title = "ERR4669752") # chr4 elevated but not on a ploidy line
SWline_closer(HSC158_bins, title = "HSC158") # chrR - right arm 1130000
SWline_closer(SRR12969478_bins, title = "SRR12969478") # +chrR - right arm 1970000-2205000
SWline_closer(SRR13587806_bins, title = "SRR13587806") # +chrR- right arm, 1615000-1790000
SWline_closer(SRR13587818_bins, title = "SRR13587818") # chr& elevated but not on a ploidy line -chrR 1905000 - end
SWline_closer(SRR13587889_bins, title = "SRR13587889") # +chr3 1410000-end
SWline_closer(SRR1554302_bins, title = "SRR1554302") # x2 +chr5 start-475000
SWline_closer(SRR21359696_bins, title = "SRR21359696") # nothing clear
SWline_closer(SRR21359704_bins, title = "SRR21359704") # +chr6
SWline_closer(SRR9073715_bins, title = "SRR9073715") #  +chr1 725000-820000
SWline_closer(SRR9073716_bins, title = "SRR9073716") # 
#SWline_closer(HSC160_bins, title = "HSC160") 
dev.off()

# ERR4669752
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "ERR4669752_depth_ploidy.pdf"), width=5, height=6, font = "Times")
par(mfrow=c(3, 1))
SWline_closer(ERR4669752_bins, title = "ERR4669752", base = 2, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(ERR4669752_bins, title = "ERR4669752", base = 3, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(ERR4669752_bins, title = "ERR4669752", base = 4, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
dev.off()

#HSC158
SWline_closer(HSC158_bins, title = "HSC158") # +chrRR
medDepth <- lapply(HSC158_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC158_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(HSC158_bins[[8]]/medDepth[[3]]*2, ylim=c(0, 5), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC158_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=226, col="red") # approx 1130000
dev.off()

#SRR12969478
medDepth <- lapply(SRR12969478_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR12969478_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR12969478_bins[[8]]/medDepth[[3]]*2, ylim=c(1, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR12969478_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=394, col="red") # approx 1970000
abline(v=441, col="red") # approx 2205000
dev.off()

#SRR13587806
medDepth <- lapply(SRR13587806_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587806_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587806_bins[[8]]/medDepth[[3]]*2, ylim=c(1, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587806_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=323, col="red") # approx 1615000
abline(v=358, col="red") # approx 1790000
dev.off()

#SRR13587818
medDepth <- lapply(SRR13587818_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587818_depth_ploidy.pdf"), width=5, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1))
SWline_closer(SRR13587818_bins, title = "SRR13587818", base = 2, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR13587818_bins, title = "SRR13587818", base = 3, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR13587818_bins, title = "SRR13587818", base = 4, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
dev.off()

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587818_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587818_bins[[8]]/medDepth[[3]]*2, ylim=c(1, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587818_bins[[8]]/medDepth[[2]]*2 < 1.5)
abline(v=381, col="red") # approx 1905000
dev.off()

#SRR13587889
medDepth <- lapply(SRR13587889_bins, function(x) quantile(x, 0.5)) 
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587889_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587889_bins[[3]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587889_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 1410000
dev.off()

#SRR1554302
medDepth <- lapply(SRR1554302_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR1554302_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR1554302_bins[[5]]/medDepth[[2]]*2, ylim=c(1, 8), main = "chr5")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR1554302_bins[[5]]/medDepth[[2]]*2 > 2.5)
abline(v=95, col="red") # approx 475000
dev.off()

#SRR21359696
medDepth <- lapply(SRR21359696_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR21359696_depth_ploidy.pdf"), width=5, height=6, font = "Times")
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1))
SWline_closer(SRR21359696_bins, title = "SRR21359696", base = 2, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR21359696_bins, title = "SRR21359696", base = 3, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR21359696_bins, title = "SRR21359696", base = 4, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
dev.off()

#SRR9073715
medDepth <- lapply(SRR9073715_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR9073715_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR9073715_bins[[1]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR9073715_bins[[1]]/medDepth[[2]]*2 > 2.5)
abline(v=145, col="red") # approx 725000
abline(v=164, col="red") # approx 820000
dev.off()

#SRR9073716
medDepth <- lapply(SRR9073716_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR9073716_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR9073716_bins[[1]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr1")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR9073716_bins[[1]]/medDepth[[2]]*2 > 2.5)
abline(v=145, col="red") # approx 725000
abline(v=164, col="red") # approx 820000
dev.off()

###############################################################################
# read in look closerlook2 data
###############################################################################
final_col_data <- cbind(initcol[,1:2], finalcol[,"SRR13587712"]) #+chrR 1795000-end
final_col_data <- as.data.frame(final_col_data)

SRR13587712_ddn <- split(final_col_data, final_col_data["Chr"])
SRR13587712_bins <- calculate_bins(SRR13587712_ddn, 3)

#SRR13587712
medDepth <- lapply(SRR13587712_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR9073716_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587712_bins[[8]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587712_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=359, col="red") # approx 1795000
dev.off()

###############################################################################
# read in hopefully the final data set
###############################################################################
#final_final_col_data <- cbind(initcol[,1:2], finalfinalcol)
final_col_data <- as.data.frame(finalfinalcol)

isol_final <- names(final_col_data)[3:length(final_col_data)]
for(i in 1:length(isol_final)){
  splitname <- paste0(isol_final[i], "_ddn")
  dataset <- final_col_data[, isol_final[i]]
  dataset <- cbind(final_col_data[,1:2], dataset)
  assign(splitname, split(dataset, dataset["Chr"]))
  covlist <- eval(parse(text = splitname))
  binname <- paste0(isol_final[i], "_bins")
  assign(binname, calculate_bins(covlist, 3))
  print(binname)
}


pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "FinalCloserIsol.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(6, 2), mar=c(1, 1, 1, 1))
SWline_closer(HSC015_bins, title = "HSC015") # +chr3 675000-785000
SWline_closer(HSC124_bins, title = "HSC124") # +chr2 2005000-2050000 (too small)
SWline_closer(HSC167_bins, title = "HSC167") # chr4 150000-205000
SWline_closer(SRR13587587_bins, title = "SRR13587587") # +chr6 2005000-2050000 (too small)
SWline_closer(SRR13587609_bins, title = "SRR13587609") # +chr6 2005000-2050000 (too small)
SWline_closer(SRR13587650_bins, title = "SRR13587650") # +chrR 1575000-1620000 (too small)
SWline_closer(SRR13587668_bins, title = "SRR13587668") # +chr6 2005000-2050000 (too small)
SWline_closer(SRR13587758_bins, title = "SRR13587758") # +chr6 2005000-2050000 (too small)
SWline_closer(SRR13587672_bins, title = "SRR13587672") # +chr6 2005000-2050000 (too small)
SWline_closer(SRR13587746_bins, title = "SRR13587746") # 4N +chr7
SWline_closer(SRR13587634_bins, title = "SRR13587634") # -chr7 770000-820000
SWline_closer(SRR13587873_bins, title = "SRR13587873") # CNV on chr3 at ALS7 - too small
SWline_closer(SRR13587879_bins, title = "SRR13587879") # +chr6 795000-970000
SWline_closer(SRR13587883_bins, title = "SRR13587883") # +chr6 795000-975000
SWline_closer(SRR1544598_bins, title = "SRR1544598") # +chr4 +chr5
SWline_closer(SRR1544616_bins, title = "SRR1544616") # +chr5
SWline_closer(SRR1548200_bins, title = "SRR1548200") # +chr3
SWline_closer(SRR1548207_bins, title = "SRR1548207") # +chr6
SWline_closer(SRR1548320_bins, title = "SRR1548320") # - 2x chrR 1905000-2285000
SWline_closer(SRR6669855_bins, title = "SRR6669855") # -chrR 485000-535000
SWline_closer(SRR6669870_bins, title = "SRR6669870") # -chrR 1-200000
SWline_closer(SRR6669886_bins, title = "SRR6669886") # chr3 1150000-1055000
SWline_closer(SRR6669890_bins, title = "SRR6669890") # -chrR 5000-245000
dev.off()

# SRR13587746
pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587746_depth_ploidy.pdf"), width=5, height=6, font = "Times")
par(mfrow=c(3, 1))
SWline_closer(SRR13587746_bins, title = "SRR13587746", base = 2, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR13587746_bins, title = "SRR13587746", base = 3, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
SWline_closer(SRR13587746_bins, title = "SRR13587746", base = 4, maxP = 6, plotYaxis = TRUE, standFormula = "chr2") # 
dev.off()

#HSC015
medDepth <- lapply(HSC015_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC015_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(HSC015_bins[[3]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC015_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=135, col="red") # approx 675000
abline(v=157, col="red") # approx 785000
dev.off()

#HSC124
medDepth <- lapply(HSC124_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC124_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(HSC124_bins[[2]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr2")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC124_bins[[2]]/medDepth[[6]]*2 > 2.5)
abline(v=401, col="red") # approx 2005000
abline(v=410, col="red") # approx 2050000
dev.off()

#SRR13587758
medDepth <- lapply(SRR13587758_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587758_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587758_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587758_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()

###
#SRR13587587
medDepth <- lapply(SRR13587587_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587587_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR13587587_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587587_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()

#SRR13587609
medDepth <- lapply(SRR13587609_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587609_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR13587609_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587609_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()

#SRR13587668
medDepth <- lapply(SRR13587668_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587668_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR13587668_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587668_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()

#SRR13587672
medDepth <- lapply(SRR13587672_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587672_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR13587672_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587672_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()

#SRR13587873
medDepth <- lapply(SRR13587873_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587873_depth.pdf"), width=4, height=6, font = "Times")
plot(SRR13587873_bins[[3]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587873_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=282, col="red") # approx 2005000
abline(v=292, col="red") # approx 2050000 - too small
dev.off()


#HSC167
medDepth <- lapply(HSC167_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "HSC167_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(HSC167_bins[[4]]/medDepth[[6]]*2, ylim=c(1, 4), main = "chr4")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(HSC167_bins[[4]]/medDepth[[2]]*2 > 2.5)
abline(v=30, col="red") # approx 150000
abline(v=41, col="red") # approx 205000 - too small
dev.off()

#SRR6669886
medDepth <- lapply(SRR6669886_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669886_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR6669886_bins[[3]]/medDepth[[2]]*2, ylim=c(1, 5), main = "chr3")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669886_bins[[3]]/medDepth[[2]]*2 > 2.5)
abline(v=211, col="red") # approx 1055000
abline(v=230, col="red") # approx 1150000 
dev.off()

#SRR13587634
medDepth <- lapply(SRR13587634_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587634_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587634_bins[[7]]/medDepth[[2]]*2, ylim=c(0.5, 4), main = "chr7")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587634_bins[[7]]/medDepth[[2]]*2 < 1.5)
abline(v=154, col="red") # approx 770000
abline(v=164, col="red") # approx 820000 
dev.off()

#SRR13587650
medDepth <- lapply(SRR13587650_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587650_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587650_bins[[8]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587650_bins[[8]]/medDepth[[2]]*2 > 2.5)
abline(v=315, col="red") # approx 1575000
abline(v=324, col="red") # approx 1620000 
dev.off()

#SRR13587879
medDepth <- lapply(SRR13587879_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587634_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587879_bins[[6]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587879_bins[[6]]/medDepth[[2]]*2 > 2.5)
abline(v=159, col="red") # approx 795000
abline(v=194, col="red") # approx 970000 
dev.off()

#SRR13587883
medDepth <- lapply(SRR13587883_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR13587883_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR13587883_bins[[6]]/medDepth[[2]]*2, ylim=c(1, 4), main = "chr6")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR13587883_bins[[6]]/medDepth[[2]]*2 > 2.5)
abline(v=159, col="red") # approx 795000
abline(v=195, col="red") # approx 975000 
dev.off()

#SRR1548320
medDepth <- lapply(SRR1548320_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR1548320_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR1548320_bins[[8]]/medDepth[[2]]*2, ylim=c(0, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR1548320_bins[[8]]/medDepth[[2]]*2 < 1.5)
abline(v=381, col="red") # approx 1905000
abline(v=457, col="red") # approx 2285000 
dev.off()

#SRR6669855
medDepth <- lapply(SRR6669855_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
        "SRR6669855_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR6669855_bins[[8]]/medDepth[[2]]*2, ylim=c(0, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669855_bins[[8]]/medDepth[[2]]*2 < 1.5)
abline(v=97, col="red") # approx 485000
abline(v=107, col="red") # approx 535000 
dev.off()

#SRR6669870
medDepth <- lapply(SRR6669870_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669870_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR6669870_bins[[8]]/medDepth[[2]]*2, ylim=c(0, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669870_bins[[8]]/medDepth[[2]]*2 < 1.5)
abline(v=1, col="red") # approx 1
abline(v=40, col="red") # approx 200000 
dev.off()

#SRR6669890
medDepth <- lapply(SRR6669890_bins, function(x) quantile(x, 0.5)) 

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "depth_singleIsolates", 
         "SRR6669890_depth.pdf"), width=4, height=6, font = "Times")
#par(mfrow=c(2, 1))
plot(SRR6669890_bins[[8]]/medDepth[[2]]*2, ylim=c(0, 4), main = "chrR")
abline(h = 1, lty=3)
abline(h = 2, lty=3)
abline(h = 3, lty=3)
abline(h = 4, lty=3)
which(SRR6669890_bins[[8]]/medDepth[[2]]*2 < 1.5)
abline(v=2, col="red") # approx 5000
abline(v=49, col="red") # approx 245000 
dev.off()