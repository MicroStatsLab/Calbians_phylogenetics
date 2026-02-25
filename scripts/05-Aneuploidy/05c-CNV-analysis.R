###############################################################################
# required libraries
###############################################################################
library(tidyverse)
library(here)
library(readxl)

###############################################################################
# functions to find chromosome lengths
###############################################################################
initcol <-  fread(here("data_in", "05-Aneuploidy","depth_1_30_coverage.txt"), sep="\t", header=T)
minData <- initcol[,1:2]
minData_ddn <- split(minData, minData$Chr)
lapply(minData_ddn, function(x) nrow(x))


centromere_positions <- data.table(
  CHROM = c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chrR"),  # Chromosomes
  start = c(1563038,1927255,823334,992579,468716,980040,425812,1743190 ), # Start of centromere region
  end = c(1565967, 1930214, 826482,996216 ,471745,983792,428712,1747664)    # End of centromere region
)

Geography_colour <- c("Africa"= "#d55e00", "Asia" = "#359e73", "Europe" = "#f748a5", "North America"="#2271b2", "South America"="#f0e442")


###############################################################################
# read CNV dataset
###############################################################################
CNV <- read_excel("data_out/05-Aneuploidy/250709PASS_CNV.xlsx")

#total = 41, each isolate gets 0.02439024 vertical space

#chr1 - 9
#chr1 - 3188548

CNV[complete.cases(CNV$Chr1),c("Isolate_Name", "Chr1")]

# 1 FDAARGOS_656 725000-820000                                                    
# 3 P94015       435000-850000                                                    
# 4 Bougn15      405000-460000, 1795000-1895000                                   
# 5 A1605/1993   1-1695000                                                        
# 6 M1653/93     480000-3025000                                                   
# 7 HSC122       460000-3188548     

#7.74x7.75
#pdf(here("Calb_manu_figures_tables", "05-Aneuploidy", "250722CNVchanges.pdf"), width = 5, height = 8)
#par(mfrow = c(6, 1), mar=c(2, 7, 1, 1), oma=c(1, 1, 1, 1))
par(mar=c(1, 5, 0, 1), oma=c(2, 1, 1, 1))
par(fig = c(0, 10, 8.451, 10)/10)
plot(c(725000, 820000), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,9), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
abline(v= 2500000, col = grey(0.85))
abline(v= 3000000, col = grey(0.85))
points(c(1, 3188548), c(7, 7), type="l")
points(c(1, 3188548), c(8, 8), type="l")
points(c(1, 3188548), c(9, 9), type="l")
points(c(725000, 820000), c(1, 1), type="l")
points(c(435000, 820000), c(4, 4), type="l")
points(c(405000, 460000), c(5, 5), type="l")
points(c(1795000, 1895000), c(5, 5), type="l")
points(c(1, 1695000), c(6, 6), type="l")
points(c(480000, 3025000), c(2, 2), type="l")
points(c(460000, 3188548), c(3, 3), type="l")
axis(2, las=2, c("8840", "K138", "LS42", "FDAARGOS_656", "P94015", "Bougn15", "A1605/1993", "M1653/93", "HSC122"), at=c(7, 8, 9, 1, 4, 5, 6, 2, 3), cex.axis=0.7, pos=0)
#axis(1, cex.axis=0.8, labels=FALSE)
arrows(0, 0.2, 3188548, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1563038, 1565967), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)
#text(c(500000, 1000000), 0, labels=c("500000", "1000000"))
mtext("Chr1", side=3, adj=0.05)

#chr2 - 2232035 - 3
CNV[complete.cases(CNV$Chr2),c("Isolate_Name", "Chr2")]
#5106         475000-2232035

par(mar=c(1, 5, 0, 1), oma=c(2, 1, 1, 1))
par(fig = c(0, 10, 8.451, 10)/10)
plot(c(1, 2232035), c(3, 3), type="l", xlim=c(0, 3400000), ylim=c(0,3), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
points(c(1, 2232035), c(3, 3), type="l")
points(c(1, 2232035), c(2, 2), type="l")
points(c(475000, 2232035), c(1, 1), type="l")
axis(2, las=2, c("Bougn15", "HAME38", "5106"), at=c(3, 2, 1), cex.axis=0.7, pos=0)
arrows(0, 0.2, 2232035, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1927255, 1930214), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)

#chr3 - 1799407 - 2
CNV[complete.cases(CNV$Chr3),c("Isolate_Name", "Chr3")]
#GZ239        1410000-1799407   
par(mar=c(1, 5, 0, 1), oma=c(2, 1, 1, 1))
par(fig = c(0, 10, 8.451, 10)/10)
plot(c(1, 1799407), c(2, 2), type="l", xlim=c(0, 3400000), ylim=c(0,3), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
points(c(1, 1799407), c(2, 2), type="l")
points(c(1410000, 1799407), c(1, 1), type="l")
axis(2, las=2, c("Z2016", "GZ239"), at=c(2, 1), cex.axis=0.7, pos=0)
arrows(0, 0.2, 2232035, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1927255, 1930214), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)



#chr4 - 1603443

CNV[complete.cases(CNV$Chr4),c("Isolate_Name", "Chr4")]
# 1 A1605/1993   1195000-1603443                  
# 2 M1653/93     410000-850000 x2, 850000 -1330000
# 3 SY3          1185000-1603443     

par(fig = c(0, 10, 7.44, 8.451)/10)
par(new=T)
plot(c(1195000, 1603443), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,3), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
points(c(1195000, 1603443), c(1, 1), type="l")
points(c(1195000, 1603443), c(1, 1), type="l")
points(c(410000, 850000), c(3, 3), type="l", lwd=2)
points(c(850000, 1330000), c(3, 3), type="l")
points(c(1185000, 1603443), c(2, 2), type="l")
axis(2, las=2, c("A1605/1993", "M1653/93", "SY3"), at=c(3, 1, 2), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1603443, 0.2, length=0, col="maroon", lwd=1.5)
points(c(992579, 996216), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr4", side=3, adj=0.05)


#chr5
CNV[complete.cases(CNV$Chr5),c("Isolate_Name", "Chr5")]
#chr5 - 1190928

# 1 J205         1-735000      
# 2 A1605/1993   1-470000      
# 3 HSC122       1-470000      
# 4 HSC139       1-470000      
# 5 3120         475000-1190928
# 6 3107         1-475000      
# 7 5106         1-390000     

par(fig = c(0, 10, 5.738, 7.44)/10)
par(new=T)
plot(c(1, 735000), c(2, 2), type="l", xlim=c(0, 3400000), ylim=c(0,7), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
points(c(1, 735000), c(2, 2), type="l")
points(c(1, 470000), c(6, 6), type="l")
points(c(1, 470000), c(5, 5), type="l")
points(c(1, 470000), c(4, 4), type="l")
points(c(475000, 1190928), c(1, 1), type="l")
points(c(1, 475000), c(3, 3), type="l")
points(c(1, 390000), c(7, 7), type="l", lwd=2)
axis(2, las=2, c("J205", "A1605/1993", "HSC122", "HSC139", "3120", "3107","5106" ), at=c(2, 6, 5, 4, 1, 3, 7), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1190928, 0.2, length=0, col="maroon", lwd=1.5)
points(c(468716, 471745), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr5", side=3, adj=0.05)

#chr6
#chr6 - 1033530

CNV[complete.cases(CNV$Chr6),c("Isolate_Name", "Chr6")]
# 1 A20          435000-1033530
# 2 J205         1-880000      
# 3 HSC131       1-865000      
# 4 A1605/1993   1- 860000     
# 5 M1653/93     330000-875000 
# 6 A_20_3       1-880000      
# 7 HSC139       205000-875000   

par(fig = c(0, 10, 4.036, 5.738)/10)
par(new=T)
plot(c(435000, 1033530), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,7), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
points(c(435000, 1033530), c(1, 1), type="l")
points(c(1, 880000), c(5, 5), type="l")
points(c(1, 865000), c(6, 6), type="l")
points(c(1, 860000), c(7, 7), type="l")
points(c(330000, 875000), c(2, 2), type="l")
points(c(1, 880000), c(4, 4), type="l")
points(c(205000, 875000), c(3, 3), type="l")
axis(2, las=2, c("A20", "J205", "HSC131", "A1605/1993", "M1653/93", "A_20_3", "HSC139"), at=1:7, cex.axis=0.7, pos=0)
arrows(0, 0.2, 1033530, 0.2, length=0, col="maroon", lwd=1.5)
points(c(980040, 983792), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr6", side=3, adj=0.05)

#chr7
#chr7 - 949616

CNV[complete.cases(CNV$Chr7),c("Isolate_Name", "Chr7")]
# A20          585000-949616
# LJ10         1-245000     
# HSC131       230000-949616
# A1605/1993   250000-949616
# M1653/93     230000-949616
par(fig = c(0, 10, 2.759, 4.036)/10)
par(new=T)
plot(c(585000, 949616), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,5), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
points(c(585000, 949616), c(1, 1), type="l")
points(c(1, 245000), c(5, 5), type="l")
points(c(230000, 949616), c(3, 3), type="l")
points(c(250000, 949616), c(4, 4), type="l")
points(c(230000, 949616), c(2, 2), type="l")
axis(2, las=2, c("A20", "LJ10", "HSC131", "A1605/1993", "M1653/93"), at=c(1, 5, 3, 4, 2), cex.axis=0.7, pos=0)
arrows(0, 0.2, 949616, 0.2, length=0, col="maroon", lwd=1.5)
points(c(425812, 428712), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr7", side=3, adj=0.05)

#chrR
#chrR - 2286389

CNV[complete.cases(CNV$Chr8),c("Isolate_Name", "Chr8")]
# 1 HSC158       1130000-2286389                                      
# 2 P60002       745000-985000, 985000-2286389 x2                     
# 3 LJ10         1355000-2286389                                      
# 4 J2010        1615000-1790000                                      
# 5 A1605/1993   75000-1765000                                        
# 6 HS78         905000 - 2286389                                     
# 7 M1653/93     235000-2286389, 2x485000-535000                      
# 8 SF53         1795000-2286389                                      
# 9 A_20_3       1435000-2286389                                      
# 10 EGPHCRM4     1160000-2286389                                      
# 11 2018-F4-39   1362000-2286389                                      
# 12 ATCC 64124   1970000-2205000                                      
# 13 SY3          1355000-2286389   
par(fig = c(0, 10, 0, 2.759)/10)
par(new=T)
plot(c(1130000, 2286389), c(9, 9), type="l", xlim=c(0, 3400000), ylim=c(0,13), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
points(c(1130000, 2286389), c(9, 9), type="l")
points(c(1130000, 2286389), c(9, 9), type="l")
points(c(745000, 985000), c(11, 11), type="l")
points(c(985000, 2286389), c(11, 11), type="l", lwd=2)
points(c(1355000, 2286389), c(5, 5), type="l")
points(c(1615000, 1790000), c(3, 3), type="l")
points(c(75000, 1765000), c(13, 13), type="l")
points(c(905000, 2286389), c(10, 10), type="l")
points(c(235000, 2286389), c(12, 12), type="l")
points(c(485000, 535000), c(12, 12), type="l", lwd=2)
points(c(1795000, 2286389), c(2, 2), type="l")
points(c(1435000, 2286389), c(4, 4), type="l")
points(c(1160000, 2286389), c(8, 8), type="l")
points(c(1362000, 2286389), c(7, 7), type="l")
points(c(1970000, 2205000), c(1, 1), type="l")
points(c(1355000, 2286389), c(6, 6), type="l")
axis(2, las=2, c("HSC158", "P60002", "LJ10", "J2010", "A1605/1993", "HS78", "M1653/93", "SF53", "A_20_3", "EGPHCRM4", "2018-F4-39", "ATCC 64124", "SY3"), at=c(9, 11, 5, 3, 13, 10, 12, 2, 4, 8, 7, 1, 6), cex.axis=0.7, pos=0)
arrows(0, 0.2, 2286389, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1743190, 1747664), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")

mtext("ChrR", side=3, adj=0.05)
axis(1, cex.axis=0.8, at = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000), labels=c("0", "500000", "1000000", "1500000", "2000000", "2500000", "3000000"), cex.axis=1.1)
mtext("Chromosome position (bp)", side=1, outer=TRUE, line=1.5)
dev.off()
