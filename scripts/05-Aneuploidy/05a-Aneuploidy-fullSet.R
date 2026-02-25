library(here)
library(data.table)

# how long is each chromosome?
#lapply(Coverage, function(x) length(t(x[,1])))

metadata <- read.csv(here("data_in","Calb_manu_tables", "TableS1.csv"), header = TRUE) #1180
metadata_pass <- subset(metadata, PASS_Het == 1) #842

# Renamed bin function
calculate_bins <- function(ddn, line, len = 5000) {
  line.ddn <- lapply(ddn, "[", c(1, 2, line))
  bin_cov <- lapply(line.ddn, function(chrom_data) {
    chrom_data$bin <- cut(chrom_data$locus, breaks = seq(0, max(chrom_data$locus), by = len), labels = FALSE)
    tapply(chrom_data[[3]], chrom_data$bin, mean, na.rm = TRUE)
  })
  names(bin_cov) <- names(ddn)
  return(bin_cov)
}

SWline <- function(bin, title){
  nm <-deparse(substitute(bin))
  # calculate the number of positions in the first chromosome
  x0 <- length(bin[[1]])
  stand <- mean(unlist(bin[[1]]))
  plot(1:x0,  2*unlist(bin[[1]])/stand,  type="l", col="red", xlim=c(0, length(unlist(bin))), xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0, 4), main = "")
  for (i in 2:8){
    x1 <- length(bin[[i]])+x0-1
    x <- x0:x1
    stand <- median(unlist(bin[[1]]))
    points(c(x0:x1), 2*unlist(bin[[i]]/stand), type="l", col=col[i])
    x0 <- x1
  }
  mtext(title, side=3, adj=0.01, cex=0.7)  
  abline(h=1, lty=2)
  abline(h=2, lty=2)
  abline(h=3, lty=2)
  abline(h=4, lty=2)
  #axis(1,at=c(49,148,254,375,509,671,862,1067,1282,1511,1761,2037,2323),cex.axis=0.5,labels=FALSE)
}

SWline_closer <- function(bin, title, base = 2, minP = 0, maxP = 4, plotYaxis = FALSE){
  nm <-deparse(substitute(bin))
  # calculate the number of positions in the first chromosome
  x0 <- length(bin[[1]])
  stand <- median(unlist(lapply(bin, function(x) median(x))))
  #stand <- mean(unlist(bin[[1]]))
  plot(1:x0,  base*unlist(bin[[1]])/stand,  type="l", col="red", xlim=c(0, length(unlist(bin))), xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(minP, maxP), main = "")
  for (i in 2:8){
    x1 <- length(bin[[i]])+x0-1
    x <- x0:x1
    stand <- median(unlist(lapply(bin, function(x) median(x))))
    points(c(x0:x1), base*unlist(bin[[i]]/stand), type="l", col=col[i])
    x0 <- x1
  }
  mtext(title, side=3, adj=0.01, cex=0.7)  
  for(k in minP:maxP) abline(h=k, lty=2)
  if (plotYaxis) axis(2, las=2)
  #axis(1,at=c(49,148,254,375,509,671,862,1067,1282,1511,1761,2037,2323),cex.axis=0.5,labels=FALSE)
}


col <- c(rep(c("red", "blue"), 7))

loaddata <- function(x) fread(here("data_in", "05-Aneuploidy", "samtools_depth", x), sep="\t", header=T)
#cov31_60 <- loaddata("depth_31_60_coverage.txt")

splitdata <- function(x){
  df <- as.data.frame(x)
  split(df, df["Chr"])
}

firstnum <- seq(1, 1111, by = 30)
secondnum <- seq(30, 1140, by = 30)
ranges <- paste(firstnum, secondnum, sep="_")

##########
#load in datasets
##########

for(i in ranges[1:9]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  splitname <- paste0("cov", i, "_ddn")
  print(splitname)
  data <- eval(parse(text=datasetname))
  assign(splitname, split(data, data["Chr"]))
  covlist <- eval(parse(text = splitname))
  isol <- names(data)[3:32]
  for(k in 1:30){
    assign(isol[k], calculate_bins(covlist, k+2))
    print(isol[k])
  }
  isol_pass <- subset(isol, isol %in% metadata_pass$Accession)
  pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", paste0(datasetname, "_coverage.pdf")), width=8, height=9, font = "Times")
  par(mfrow=c(8, ceiling(length(isol_pass)/8)), mar=c(1, 1, 1, 1))
  for(j in 1:length(isol_pass)){
    SWline(eval(parse(text = isol_pass[j])), title = isol_pass[j])   
  }
  dev.off()
}

for(i in ranges[10:18]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  splitname <- paste0("cov", i, "_ddn")
  print(splitname)
  data <- eval(parse(text=datasetname))
  assign(splitname, split(data, data["Chr"]))
  covlist <- eval(parse(text = splitname))
  isol <- names(data)[3:32]
  for(k in 1:30){
    assign(isol[k], calculate_bins(covlist, k+2))
    print(isol[k])
  }
  isol_pass <- subset(isol, isol %in% metadata_pass$Accession)
  pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", paste0(datasetname, "_coverage.pdf")), width=8, height=9, font = "Times")
  par(mfrow=c(8, ceiling(length(isol_pass)/8)), mar=c(1, 1, 1, 1))
  for(j in 1:length(isol_pass)){
    SWline(eval(parse(text = isol_pass[j])), title = isol_pass[j])   
  }
  dev.off()
}

for(i in ranges[19:27]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  splitname <- paste0("cov", i, "_ddn")
  print(splitname)
  data <- eval(parse(text=datasetname))
  assign(splitname, split(data, data["Chr"]))
  covlist <- eval(parse(text = splitname))
  isol <- names(data)[3:32]
  for(k in 1:30){
    assign(isol[k], calculate_bins(covlist, k+2))
    print(isol[k])
  }
  isol_pass <- subset(isol, isol %in% metadata_pass$Accession)
  if(length(isol_pass) > 0){
    pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", paste0(datasetname, "_coverage.pdf")), width=8, height=9, font = "Times")
    par(mfrow=c(8, ceiling(length(isol_pass)/8)), mar=c(1, 1, 1, 1))
    for(j in 1:length(isol_pass)){
      SWline(eval(parse(text = isol_pass[j])), title = isol_pass[j])   
    }
    dev.off()
  }
  rm(datasetname)
  rm(splitname)
}

for(i in ranges[28:36]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  splitname <- paste0("cov", i, "_ddn")
  print(splitname)
  data <- eval(parse(text=datasetname))
  assign(splitname, split(data, data["Chr"]))
  covlist <- eval(parse(text = splitname))
  isol <- names(data)[3:32]
  for(k in 1:30){
    assign(isol[k], calculate_bins(covlist, k+2))
    print(isol[k])
  }
  isol_pass <- subset(isol, isol %in% metadata_pass$Accession)
  if(length(isol_pass) > 0){
    pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", paste0(datasetname, "_coverage.pdf")), width=8, height=9, font = "Times")
    par(mfrow=c(8, ceiling(length(isol_pass)/8)), mar=c(1, 1, 1, 1))
    for(j in 1:length(isol_pass)){
      SWline(eval(parse(text = isol_pass[j])), title = isol_pass[j])   
    }
    dev.off()
  }
  rm(datasetname)
  rm(splitname)
}

for(i in ranges[37:38]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  splitname <- paste0("cov", i, "_ddn")
  print(splitname)
  data <- eval(parse(text=datasetname))
  assign(splitname, split(data, data["Chr"]))
  covlist <- eval(parse(text = splitname))
  isol <- names(data)[3:length(data)]
  for(k in 1:length(isol)){
    assign(isol[k], calculate_bins(covlist, k+2))
    print(isol[k])
  }
  isol_pass <- subset(isol, isol %in% metadata_pass$Accession)
  if(length(isol_pass) > 0){
    pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", paste0(datasetname, "_coverage.pdf")), width=8, height=9, font = "Times")
    par(mfrow=c(8, ceiling(length(isol_pass)/8)), mar=c(1, 1, 1, 1))
    for(j in 1:length(isol_pass)){
      SWline(eval(parse(text = isol_pass[j])), title = isol_pass[j])   
    }
    dev.off()
  }
  rm(datasetname)
  rm(splitname)
}

#closer look - scored manually from graphs produced above
lookclose <- c("ERR2708450", "ERR2708456", "ERR4669765", "ERR4669777", "HSC046", "HSC050", "HSC083", "HSC122", "HSC131", "HSC139", "HSC158", "HSC160", "HSC182", "HSC210", "SRR10209617", "SRR13587535", "SRR13587542", "SRR13587569", "SRR13587577", "SRR13587579", "SRR13587596", "SRR13587612", "SRR13587639", "SRR13587643", "SRR13587657", "SRR13587673", "SRR13587680", "SRR13587748", "SRR13587782", "SRR13587761", "SRR13587868", "SRR13587878", "SRR13587889", "SRR13741117", "SRR13741118", "SRR13741121", "SRR1553997", "SRR1554177", "SRR1554178", "SRR1554287", "SRR1554289", "SRR1554296", "SRR1554297", "SRR1554304", "SRR1554310", "SRR16959048", "SRR19696103", "SRR23500659", "SRR3593469", "SRR392813", "SRR393519", "SRR397731", "SRR530262", "SRR543720", "SRR543723", "SRR543724", "SRR640897", "SRR6669863", "SRR6669867", "SRR6669868", "SRR6669891", "SRR6669909", "SRR6669920", "SRR6669922", "SRR6669940", "SRR6669942", "SRR6669948", "SRR6669989", "SRR6670027", "SRR13587749")
length(lookclose)/nrow(metadata_pass)

# pull out the closer look datasets
for(i in ranges[1:9]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  datasetname_sub <- paste0("cov", i, "_sub")
  data <- eval(parse(text=datasetname))
  data_sub <- data[, names(data)[names(data) %in% lookclose]]
  assign(datasetname_sub, data_sub)
}

lookclose_data_1_270 <- cbind(cov1_30_sub, cov31_60_sub, cov61_90_sub, cov91_120_sub, cov121_150_sub, cov151_180_sub, cov181_210_sub, cov211_240_sub, cov241_270_sub)

write.csv(lookclose_data_1_270, here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_1_270.csv"))

for(i in ranges[10:19]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  datasetname_sub <- paste0("cov", i, "_sub")
  data <- eval(parse(text=datasetname))
  data_sub <- data[, names(data)[names(data) %in% lookclose]]
  assign(datasetname_sub, data_sub)
}

lookclose_data_271_570 <- cbind(cov271_300_sub, cov331_360_sub, cov361_390_sub, cov451_480_sub, cov481_510_sub, cov541_570_sub)
names(lookclose_data_271_570)[4] <- names(cov331_360)[names(cov331_360) %in% lookclose]

write.csv(lookclose_data_271_570, here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_271_570.csv"))

for(i in ranges[20:29]){
  datasetname <- paste0("cov", i)
  filename <- paste0("depth_", i, "_coverage.txt")
  assign(datasetname, as.data.frame(loaddata(filename)))
  datasetname_sub <- paste0("cov", i, "_sub")
  data <- eval(parse(text=datasetname))
  data_sub <- data[, names(data)[names(data) %in% lookclose]]
  assign(datasetname_sub, data_sub)
}

lookclose_data_571_870 <- as.data.frame(cbind(cov571_600_sub, cov751_780_sub)) #2
names(lookclose_data_571_870)[1] <- names(cov571_600)[names(cov571_600) %in% lookclose]
names(lookclose_data_571_870)[2] <- names(cov751_780)[names(cov751_780) %in% lookclose]

write.csv(lookclose_data_571_870, here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_571_870.csv"))

for(i in ranges[30:38]){
  datasetname <- paste0("cov", i)
  print(datasetname)
  filename <- paste0("depth_", i, "_coverage.txt")
  print(filename)
  assign(datasetname, as.data.frame(loaddata(filename)))
  datasetname_sub <- paste0("cov", i, "_sub")
  data <- eval(parse(text=datasetname))
  data_sub <- data[, names(data)[names(data) %in% lookclose]]
  assign(datasetname_sub, data_sub)
}

lookclose_data_871_1140 <- cbind(cov871_900_sub, cov901_930_sub, cov931_960_sub, cov961_990_sub, cov991_1020_sub, cov1051_1080_sub, cov1081_1110_sub, cov1111_1140_sub) #22
names(lookclose_data_871_1140)[20] <- names(cov1051_1080)[names(cov1051_1080) %in% lookclose]
names(lookclose_data_871_1140)[21] <- names(cov1081_1110)[names(cov1081_1110) %in% lookclose]
names(lookclose_data_871_1140)[22] <- names(cov1111_1140)[names(cov1111_1140) %in% lookclose]

write.csv(lookclose_data_871_1140, here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_871_1140.csv"))

#read in look closer data
initcol <-  fread(here("data_in", "05-Aneuploidy", "samtools_depth", "depth_1_30_coverage.txt"), sep="\t", header=T)
lookclose_data_1_270 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_1_270.csv"), sep=",", header=T)  
lookclose_data_1_270 <- cbind(initcol[,1:2], lookclose_data_1_270)
lookclose_data_1_270 <- as.data.frame(lookclose_data_1_270)

#
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

lookclose_data_271_570 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_271_570.csv"), sep=",", header=T)  
lookclose_data_271_270 <- cbind(initcol[,1:2], lookclose_data_271_570)
lookclose_data_271_270 <- as.data.frame(lookclose_data_271_270)

#
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

lookclose_data_571_870 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_571_870.csv"), sep=",", header=T)  
lookclose_data_571_870 <- cbind(initcol[,1:2], lookclose_data_571_870)
lookclose_data_571_870 <- as.data.frame(lookclose_data_571_870)

#
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

lookclose_data_871_1140 <- fread(here("data_out", "05-Aneuploidy", "PASS_aneuploids", "lookclose_data_871_1140.csv"), sep=",", header=T)  
lookclose_data_871_1140 <- cbind(initcol[,1:2], lookclose_data_871_1140)
lookclose_data_871_1140 <- as.data.frame(lookclose_data_871_1140)

#
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


pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook1_20.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
SWline_closer(ERR2708450_bins, title = "ERR2708450")
SWline_closer(ERR2708456_bins, title = "ERR2708456")
SWline_closer(ERR4669765_bins, title = "ERR4669765")
SWline_closer(ERR4669777_bins, title = "ERR4669777")
SWline_closer(HSC046_bins, title = "HSC046")
SWline_closer(HSC050_bins, title = "HSC050")
SWline_closer(HSC083_bins, title = "HSC083")
SWline_closer(HSC122_bins, title = "HSC122")
SWline_closer(HSC131_bins, title = "HSC131")
SWline_closer(HSC131_bins, title = "HSC139")

SWline_closer(HSC158_bins, title = "HSC158")
SWline_closer(HSC160_bins, title = "HSC160")
SWline_closer(HSC182_bins, title = "HSC182")
SWline_closer(HSC210_bins, title = "HSC210")
SWline_closer(SRR10209617_bins, title = "SRR10209617")
SWline_closer(SRR13587535_bins, title = "SRR13587535")
SWline_closer(SRR13587542_bins, title = "SRR13587542")
SWline_closer(SRR13587569_bins, title = "SRR13587569")
SWline_closer(SRR13587577_bins, title = "SRR13587577")
SWline_closer(SRR13587579_bins, title = "SRR13587579")
dev.off()

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "HSC160_varPloidy.pdf"), width=4, height=8, font = "Times")
par(mfrow=c(6, 1), mar=c(1, 3, 1, 1))
SWline_closer(HSC160_bins, "HSC160, base = 1", base = 1, 0, 4, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 2", base = 2, 0, 4, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 3", base = 3, 0, 6, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 4", base = 4, 2, 8, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 5", base = 5, 4, 10, plotYaxis = TRUE)
SWline_closer(HSC160_bins, "HSC160, base = 6", base = 6, 4, 10, plotYaxis = TRUE)
dev.off()

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook21_40.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
SWline_closer(ERR2708450_bins, title = "SRR1554289")
SWline_closer(ERR2708456_bins, title = "SRR1554296")
SWline_closer(ERR4669765_bins, title = "SRR1554297")
SWline_closer(ERR4669777_bins, title = "SRR1554304")
SWline_closer(HSC046_bins, title = "SRR1554310")
SWline_closer(HSC050_bins, title = "SRR16959048")
#isol 571
SWline_closer(HSC083_bins, title = "SRR13587680")
SWline_closer(HSC122_bins, title = "SRR13587748")
SWline_closer(HSC131_bins, title = "SRR13587761")
SWline_closer(HSC131_bins, title = "SRR13587782")

SWline_closer(HSC158_bins, title = "SRR13587868")
SWline_closer(HSC160_bins, title = "SRR13587878")
SWline_closer(HSC182_bins, title = "SRR13587889")
SWline_closer(HSC210_bins, title = "SRR13741117")
SWline_closer(SRR10209617_bins, title = "SRR13741118")
SWline_closer(SRR13587535_bins, title = "SRR13741121")
SWline_closer(SRR13587542_bins, title = "SRR1553997")
SWline_closer(SRR13587569_bins, title = "SRR1554177")
SWline_closer(SRR13587577_bins, title = "SRR1554178")
SWline_closer(SRR13587579_bins, title = "SRR1554287")
dev.off()


pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook41_60.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
#isol271-570
SWline_closer(SRR1554289_bins, title = "SRR1554289")
SWline_closer(SRR1554296_bins, title = "SRR1554296")
SWline_closer(SRR1554297_bins, title = "SRR1554297")
SWline_closer(SRR1554304_bins, title = "SRR1554304")
SWline_closer(SRR1554310_bins, title = "SRR1554310")
SWline_closer(SRR16959048_bins, title = "SRR16959048")
#isol571-870
SWline_closer(SRR19696103_bins, title = "SRR19696103")
SWline_closer(SRR23500659_bins, title = "SRR23500659")
#isol870-
SWline_closer(SRR3593469_bins, title = "SRR3593469")
SWline_closer(SRR392813_bins, title = "SRR392813")
SWline_closer(SRR393519_bins, title = "SRR393519")
SWline_closer(SRR397731_bins, title = "SRR397731")
SWline_closer(SRR530262_bins, title = "SRR530262")
SWline_closer(SRR543720_bins, title = "SRR543720")
SWline_closer(SRR543723_bins, title = "SRR543723")
SWline_closer(SRR543724_bins, title = "SRR543724")
SWline_closer(SRR640897_bins, title = "SRR640897")
SWline_closer(SRR6669863_bins, title = "SRR6669863")
SWline_closer(SRR6669867_bins, title = "SRR6669867")
SWline_closer(SRR6669868_bins, title = "SRR6669868")
dev.off()

pdf(here("Calb_manu_figures_tables", "05-aneuploidy", "samtools_depth", "CloserLook61_70.pdf"), width=8, height=11, font = "Times")
par(mfrow=c(10, 2), mar=c(1, 1, 1, 1))
SWline_closer(SRR6669891_bins, title = "SRR6669891")
SWline_closer(SRR6669909_bins, title = "SRR6669909")
SWline_closer(SRR6669920_bins, title = "SRR6669920")
SWline_closer(SRR6669922_bins, title = "SRR6669922")
SWline_closer(SRR6669940_bins, title = "SRR6669940")
SWline_closer(SRR6669942_bins, title = "SRR6669942")
SWline_closer(SRR6669948_bins, title = "SRR6669948")
SWline_closer(SRR6669989_bins, title = "SRR6669989")
SWline_closer(SRR6670027_bins, title = "SRR6670027")
SWline_closer(SRR13587749_bins, title = "SRR13587749")
dev.off()


