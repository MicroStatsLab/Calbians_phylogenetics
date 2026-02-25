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

# centromere
centromere_positions <- data.table(
  CHROM = c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chrR"),  # Chromosomes
  start = c(1563038,1927255,823334,992579,468716,980040,425812,1743190 ), # Start of centromere region
  end = c(1565967, 1930214, 826482,996216 ,471745,983792,428712,1747664)    # End of centromere region
)

Geography_colour <- c("Africa"= "#d55e00", "Asia" = "#359e73", "Europe" = "#f748a5", "North America"="#2271b2", "South America"="#f0e442")


###############################################################################
# read datasets
###############################################################################
CNV <- read_excel("data_out/05-Aneuploidy/250808PASS_CNV.xlsx")
Aneup <- read.csv("data_out/05-Aneuploidy/250808PASS_aneuploidy.csv")


# How many CNV/Aneup for each chromosome?
CNV[complete.cases(CNV$Chr1),c("Isolate_Name", "Chr1")] #6 (I thinkk the last two are the same)
# 1 HSC122       460000-3188548                
# 2 P94015       435000-850000                 
# 3 M1653/93     480000-3025000                
# 4 A1605/1993   1-1195000                     
# 5 Bougn15      405000-460000, 1795000-1895000
# 6 FDAARGOS_656 725000-820000                 
# 7 FDAARGOS_656 725000-820000   
Aneup[complete.cases(Aneup$Chr1),c("Isolate_Name", "Chr1")] #3
# 34          LS42    1
# 48          K138    1
# 101         8840    1
#chr1 total: 9

CNV[complete.cases(CNV$Chr2),c("Isolate_Name", "Chr2")] #1
# 1 5106         2x'475000-2232035
Aneup[complete.cases(Aneup$Chr2),c("Isolate_Name", "Chr2")] #6
# 19       HAME38    1
# 32         QD60    1
# 34         LS42    1
# 38        Z2016    1
# 88     M1653/93    1
# 95      Bougn15    1
#chr2 total: 7

CNV[complete.cases(CNV$Chr3),c("Isolate_Name", "Chr3")] #3
#   1 HSC015       675000-785000  
# 2 GZ239        1410000-1799407
# 3 A1622/1993   1150000-1055000
Aneup[complete.cases(Aneup$Chr3),c("Isolate_Name", "Chr3")] #2
# 38        Z2016    1
# 59        TWTC6    1
#chr3 total: 5

CNV[complete.cases(CNV$Chr4),c("Isolate_Name", "Chr4")] #4
# 1 HSC167       150000-205000                    
# 2 SY3          1185000-1603443                  
# 3 M1653/93     410000-850000 x2, 850000 -1330000
# 4 A1605/1993   1210000-1603443        
Aneup[complete.cases(Aneup$Chr4),c("Isolate_Name", "Chr4")] #12
# 3    2018-F1-60    1
# 9        HSC083    1
# 41         GX22    1
# 42         GX16    1
# 55       IDB104    1
# 57         3731    1
# 62       TWTC13   -1
# 65         1649    1
# 66         2501    1
# 71         5106    1
# 78          12C    1
# 84       P78042    1
#chr4 total: 16

CNV[complete.cases(CNV$Chr5),c("Isolate_Name", "Chr5")] #7
# 1 HSC122       1-470000      
# 2 HSC139       1-470000      
# 3 J205         -1-735000     
# 4 3120         1-475000      
# 5 3281         475000-1190928
# 6 5106         2x1-390000    
# 7 A1605/1993   1-470000      
Aneup[complete.cases(Aneup$Chr5),c("Isolate_Name", "Chr5")] #15
# 48         K138    1
# 50        GZ243    1
# 57         3731    1
# 58         3733    1
# 63         1619    1
# 64         2823    1
# 65         1649    1
# 67         3034    1
# 68         3107    1
# 79       P75010    1
# 81       P60002    1
# 88     M1653/93    1
# 93      CEC4495    1
# 97           M7    1
# 98           M6    1
#chr5 total: 22

CNV[complete.cases(CNV$Chr6),c("Isolate_Name", "Chr6")] #9
# 1 HSC131       1-865000        
# 2 HSC139       3x 205000-875000
# 3 J205         -1-880000       
# 4 H11          795000-970000   
# 5 GY015        795000-970000   
# 6 A_20_3       -1-880000       
# 7 A20          435000-1033530  
# 8 M1653/93     330000-875000   
# 9 A1605/1993   1- 860000     
Aneup[complete.cases(Aneup$Chr6),c("Isolate_Name", "Chr6")] #10
# 15        HSC160    1
# 25         G1101    1
# 38         Z2016    1
# 48          K138    1
# 60         TWTC8    1
# 62        TWTC13   -1
# 75         831-2    1
# 81        P60002    1
# 89    A1605/1993    1
# 100       Bougn9    1
#chr6 total: 19

CNV[complete.cases(CNV$Chr7),c("Isolate_Name", "Chr7")] #6
# 1 HSC131       230000-949616 
# 2 LJ10         1-245000      
# 3 Z1106        -770000-820000
# 4 A20          585000-949616 
# 5 M1653/93     230000-949616 
# 6 A1605/1993   250000-949616 
Aneup[complete.cases(Aneup$Chr7),c("Isolate_Name", "Chr7")] #23
# 4   2018-F2-279    1
# 7        HSC046    1
# 8        HSC050    1
# 9        HSC083    1
# 13       HSC139    2
# 15       HSC160    1
# 17       HSC182    1
# 18       HSC210   -1
# 21         SD79    1
# 25        G1101    1
# 27         Y208    1
# 29        FJ178    1
# 34         LS42    1
# 37         Z419    1
# 44        GH257    1
# 45         W322    1
# 48         K138    2
# 50        GZ243    1
# 62       TWTC13   -1
# 80          19F    1
# 82          L26    1
# 94      CEC4484    1
# 99          MAN    1
#chr7 total: 29

CNV[complete.cases(CNV$Chr8),c("Isolate_Name", "Chr8")] #18
# 1 2018-F4-39   1362000-2286389                 
# 2 HSC158       1130000-2286389                 
# 3 ATCC 64124   1970000-2205000                 
# 4 LJ10         1355000-2286389                 
# 5 SY3          1355000-2286389                 
# 6 LS23         1575000-1620000                 
# 7 SF53         1795000-2286389                 
# 8 J2010        1615000-1790000                 
# 9 HS78         -1905000 - end                  
# 10 TWTC11       -2x 1905000-2285000             
# 11 A_20_3       1435000-2286389                 
# 12 P60002       745000-985000, 985000-2286389 x2
# 13 H4           -485000-535000                  
# 14 EGPHCRM4     1160000-2286389                 
# 15 M1653/93     235000-2286389, 2x485000-535000 
# 16 A1605/1993   75000-1765000                   
# 17 M8684/93     1-200000                        
# 18 CEC4854      - 5000-245000  
Aneup[complete.cases(Aneup$Chr8),c("Isolate_Name", "Chr8")] #5
# 23         Y265    1
# 38        Z2016    1
# 55       IDB104    1
# 93      CEC4495    1
# 96          C78    1
#chr8 total: 23

karyoSumTable <- data.frame(Chr = paste0("Chr", c(1:7, "R")), Aneuploidy = c(3, 6, 2, 12, 15, 10, 23, 5), CNV = c(6, 1, 3, 4, 7, 9, 6, 18))
karyoSumTable$Total <- karyoSumTable$Aneuploidy + karyoSumTable$CNV
plot(karyoSumTable$Aneuploidy, karyoSumTable$CNV)
cor.test(karyoSumTable$Aneuploidy, karyoSumTable$CNV)

#7.74x7.75
#pdf(here("Calb_manu_figures_tables", "05-Aneuploidy", "250722CNVchanges.pdf"), width = 5, height = 8)
#par(mfrow = c(6, 1), mar=c(2, 7, 1, 1), oma=c(1, 1, 1, 1))

#total = 9+3+2+12+21+14+27+18 = 106
1/62*10 #each isolate gets 0.1612903 (= 54 CNVs + 1 for each Chr full length aneup)
#first half = 26 total
7*0.161 # = 1.127 #chr1 
chr1place <- c(0, 5, 8.873, 10)
2*0.161 #= 0.507 #chr2 # increased
chr2place <- c(0, 5, 8.00, 8.873)
4*0.161 #= 0.644 #chr3 #increased
chr3place <- c(0, 5, 6.75, 8.00)
5*0.161 #= 2.028 #chr4
chr4place <- c(0, 5, 4.722, 6.75)
8*0.161 # = 3.549 
chr5place <- c(0, 5, 1.173, 4.722)

# second half = 35 total
10*0.161 #= 2.366 #chr6
chr6place <- c(5, 10, 7.634, 10)
7*0.161 #= 4.563 #chr7
chr7place <- c(5, 10, 3.051, 7.634)
18*0.161 #= 3.042 #chrR
chrRplace <- c(5, 10, 0.01, 3.051)

#chr1 - 7
#pdf(here("Calb_manu_figures_tables", "05-Aneuploidy", "250808Karyotypechanges.pdf"), width = 8, height = 10)
par(mar=c(1, 5, 0, 1), oma=c(2, 1, 2, 1))
par(fig = chr1place/10)
plot(c(725000, 820000), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,9), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
abline(v= 2500000, col = grey(0.85))
abline(v= 3000000, col = grey(0.85))
# 6 A1605/1993   1-1195000                     
# 5 Bougn15      405000-460000, 1795000-1895000
# 4 P94015       435000-850000  
# 3 HSC122       460000-3188548                
# 2 M1653/93     480000-3025000                
# 1 FDAARGOS_656 725000-820000                 

points(c(1, 3188548), c(7, 7), type="l")
points(c(725000, 820000), c(1, 1), type="l")
points(c(480000, 3025000), c(2, 2), type="l")
points(c(460000, 3188548), c(3, 3), type="l")
points(c(435000, 820000), c(4, 4), type="l")
points(c(405000, 460000), c(5, 5), type="l")
points(c(1795000, 1895000), c(5, 5), type="l")
points(c(1, 1195000), c(6, 6), type="l")

axis(2, las=2, c("FDAARGOS_656", "M1653/93", "HSC122", "P94015", "Bougn15", "A1605/1993", "4 isolates"), at=c(1:7), cex.axis=0.7, pos=0)
#axis(1, cex.axis=0.8, labels=FALSE)
arrows(0, 0.2, 3188548, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1563038, 1565967), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)
#text(c(500000, 1000000), 0, labels=c("500000", "1000000"))
mtext("Chr1", side=3, adj=0.05, cex = 0.85)

#chr2 - 2
CNV[complete.cases(CNV$Chr2),c("Isolate_Name", "Chr2")]
#5106         475000-2232035
nrow(Aneup[complete.cases(Aneup$Chr2),c("Isolate_Name", "Chr2")]) #6

#par(fig = c(0, 5, 7.795, 8.302)/10)
par(fig = chr2place/10)
par(new=T)
plot(c(1, 2232035), c(3, 3), type="l", xlim=c(0, 3400000), ylim=c(0,2), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
points(c(1, 2232035), c(2, 2), type="l")
points(c(475000, 2232035), c(1, 1), type="l", lwd=2)
axis(2, las=2, c("5106", "6 isolates"), at=c(1:2), cex.axis=0.7, pos=0)
arrows(0, 0.2, 2232035, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1927255, 1930214), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)
mtext("Chr2", side=3, adj=0.05, cex = 0.85)

#chr3 - 1799407 - 4
nrow(Aneup[complete.cases(Aneup$Chr3),c("Isolate_Name", "Chr3")]) #2

CNV[complete.cases(CNV$Chr3),c("Isolate_Name", "Chr3")]
# 3 HSC015       675000-785000  
# 2 A1622/1993   1150000-1055000  
# 1 GZ239        1410000-1799407

par(fig = chr3place/10)
par(new=T)
plot(c(1, 1799407), c(4, 4), type="l", xlim=c(0, 3400000), ylim=c(0,4), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
points(c(1, 1799407), c(4, 4), type="l")
points(c(675000, 785000), c(3, 3), type="l")
points(c(1150000, 1055000), c(2, 2), type="l")
points(c(1410000, 1799407), c(1, 1), type="l")
axis(2, las=2, c("GZ239", "A1622/1993", "HSC015", "2 isolates"), at=c(1:4), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1799407, 0.2, length=0, col="maroon", lwd=1.5)
points(c(823334, 826482), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon",)
mtext("Chr3", side=3, adj=0.05, cex = 0.85)


#chr4 - 1603443 - 5
CNV[complete.cases(CNV$Chr4),c("Isolate_Name", "Chr4")]
# 4 HSC167       150000-205000                    
# 3 M1653/93     410000-850000 x2, 850000 -1330000
# 2 A1605/1993   1210000-1603443 
# 1 SY3          1185000-1603443                  
nrow(Aneup[complete.cases(Aneup$Chr4),c("Isolate_Name", "Chr4")])

par(fig = chr4place/10)
par(new=T)
plot(c(1185000, 1603443), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0, 5), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
points(c(1, 1603443), c(5, 5), type="l")
points(c(1185000, 1603443), c(1, 1), type="l")
points(c(1210000, 1603443), c(2, 2), type="l")
points(c(410000, 850000), c(3, 3), type="l", lwd=2)
points(c(850000, 1330000), c(3, 3), type="l")
points(c(150000, 205000), c(4, 4), type="l", lwd=2)

axis(2, las=2, c("SY3", "A1605/1993", "M1653/93", "HSC167", "12 isolates"), at=c(1:5), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1603443, 0.2, length=0, col="maroon", lwd=1.5)
points(c(992579, 996216), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr4", side=3, adj=0.05, cex=0.85)


#chr5 - 8
CNV[complete.cases(CNV$Chr5),c("Isolate_Name", "Chr5")]
# 7 5106         2x1-390000    
# 6 HSC122       1-470000      
# 5 HSC139       1-470000      
# 4 A1605/1993   1-470000   
# 3 3120         1-475000 
# 2 J205         -1-735000     
# 1 3281         475000-1190928
nrow(Aneup[complete.cases(Aneup$Chr5),c("Isolate_Name", "Chr5")]) #15


par(fig = chr5place/10)
par(new=T)
plot(c(1, 475000), c(3, 3), type="l", xlim=c(0, 3400000), ylim=c(0,8), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
points(c(1, 1190928), c(8, 8), type="l")
points(c(1, 390000), c(7, 7), type="l", lwd=2)
points(c(1, 470000), c(6, 6), type="l")
points(c(1, 470000), c(5, 5), type="l")
points(c(1, 470000), c(4, 4), type="l")
points(c(1, 475000), c(3, 3), type="l")
points(c(1, 735000), c(2, 2), type="l", lty = 3)
points(c(475000, 1190928), c(1, 1), type="l")
axis(2, las=2, c("3281", "J205", "3120","A1605/1993" ,"HSC139","HSC122","5106", "15 isolates" ), at=c(1:8), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1190928, 0.2, length=0, col="maroon", lwd=1.5)
points(c(468716, 471745), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr5", side=3, adj=0.05, cex=0.85)

#chr6 - 10

CNV[complete.cases(CNV$Chr6),c("Isolate_Name", "Chr6")]
# 1 HSC131       1-865000        
# 2 HSC139       3x 205000-875000
# 3 J205         -1-880000       
# 4 H11          795000-970000   
# 5 GY015        795000-970000   
# 6 A_20_3       -1-880000       
# 7 A20          435000-1033530  
# 8 M1653/93     330000-875000   
# 9 A1605/1993   1- 860000    

par(fig = chr6place/10)
par(new=T)
plot(c(435000, 1033530), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,14), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
points(c(1, 1033530), c(14, 14), type="l")
points(c(1, 1033530), c(13, 13), type="l")
points(c(1, 1033530), c(12, 12), type="l")
points(c(1, 1033530), c(11, 11), type="l")
points(c(1, 1033530), c(10, 10), type="l")
points(c(1, 1033530), c(9, 9), type="l")
points(c(1, 1033530), c(8, 8), type="l")
points(c(435000, 1033530), c(1, 1), type="l")
points(c(1, 880000), c(5, 5), type="l")
points(c(1, 865000), c(6, 6), type="l")
points(c(1, 860000), c(7, 7), type="l")
points(c(330000, 875000), c(2, 2), type="l")
points(c(1, 880000), c(4, 4), type="l")
points(c(205000, 875000), c(3, 3), type="l")
axis(2, las=2, c("831-2", "A1605/1993", "Bougn9", "G1101", "HSC160", "K138", "Z2016",  "A20", "J205", "HSC131", "A1605/1993", "M1653/93", "A_20_3", "HSC139"), at=c(14:8, 1, 5, 6, 7, 2, 4, 3), cex.axis=0.7, pos=0)
arrows(0, 0.2, 1033530, 0.2, length=0, col="maroon", lwd=1.5)
points(c(980040, 983792), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
mtext("Chr6", side=3, adj=0.05, cex=0.85)

#chr7 - 27
#chr7 - 949616

CNV[complete.cases(CNV$Chr7),c("Isolate_Name", "Chr7")]
# A20          585000-949616
# LJ10         1-245000     
# HSC131       230000-949616
# A1605/1993   250000-949616
# M1653/93     230000-949616
par(fig = chr7place/10)
par(new=T)
plot(c(585000, 949616), c(1, 1), type="l", xlim=c(0, 3400000), ylim=c(0,27), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
points(c(1, 949616), c(27, 27), type="l")
points(c(1, 949616), c(26, 26), type="l")
points(c(1, 949616), c(25, 25), type="l")
points(c(1, 949616), c(24, 24), type="l")
points(c(1, 949616), c(23, 23), type="l")
points(c(1, 949616), c(22, 22), type="l")
points(c(1, 949616), c(21, 21), type="l")
points(c(1, 949616), c(20, 20), type="l")
points(c(1, 949616), c(19, 19), type="l")
points(c(1, 949616), c(18, 18), type="l")
points(c(1, 949616), c(17, 17), type="l")
points(c(1, 949616), c(16, 16), type="l")
points(c(1, 949616), c(15, 15), type="l")
points(c(1, 949616), c(14, 14), type="l")
points(c(1, 949616), c(13, 13), type="l")
points(c(1, 949616), c(12, 12), type="l")
points(c(1, 949616), c(11, 11), type="l")
points(c(1, 949616), c(10, 10), type="l")
points(c(1, 949616), c(9, 9), type="l")
points(c(1, 949616), c(8, 8), type="l")
points(c(1, 949616), c(7, 7), type="l")
points(c(1, 949616), c(6, 6), type="l")
points(c(585000, 949616), c(1, 1), type="l")
points(c(1, 245000), c(5, 5), type="l")
points(c(230000, 949616), c(3, 3), type="l")
points(c(250000, 949616), c(4, 4), type="l")
points(c(230000, 949616), c(2, 2), type="l")
axis(2, las=2, c("19F", "2018-F2-279", "CEC4484", "FJ178", "G1101", "GH257", "GZ243", "HSC046", "HSC050", "HSC083", "HSC139", "HSC160", "HSC182", "K138", "L26", "LS42", "MAN", "P60002", "SD79", "W322", "Y208", "Z419", "A20", "LJ10", "HSC131", "A1605/1993", "M1653/93"), at=c(27:6, 1, 5, 3, 4, 2), cex.axis=0.7, pos=0)
arrows(0, 0.2, 949616, 0.2, length=0, col="maroon", lwd=1.5)
points(c(425812, 428712), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")
axis(1, cex.axis=0.8, at = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000), labels=c("0", "500", "1000", "1500", "2000", "2500", "3000"), cex.axis=0.8)
mtext("Chr7", side=3, adj=0.05, cex=0.85)

#chrR
#chrR - 2286389 - 18

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
#TWTC11           - 2x 1905000-2285000
par(fig = chrRplace/10)
par(new=T)
plot(c(1130000, 2286389), c(9, 9), type="l", xlim=c(0, 3400000), ylim=c(0,18), yaxt= "n", ylab="", xaxt = "n",  bty="n")
abline(v= 500000, col = grey(0.85))
abline(v= 1000000, col = grey(0.85))
abline(v= 1500000, col = grey(0.85))
abline(v= 2000000, col = grey(0.85))
points(c(1, 2286389), c(18, 18), type="l")
points(c(1, 2286389), c(17, 17), type="l")
points(c(1, 2286389), c(16, 16), type="l")
points(c(1, 2286389), c(15, 15), type="l")
points(c(1, 2286389), c(14, 14), type="l")
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
axis(2, las=2, c("C78", "CEC4495", "IDB104", "Y265", "Z2016", "HSC158", "P60002", "LJ10", "J2010", "A1605/1993", "HS78", "M1653/93", "SF53", "A_20_3", "EGPHCRM4", "2018-F4-39", "ATCC 64124", "SY3"), at=c(18:14, 9, 11, 5, 3, 13, 10, 12, 2, 4, 8, 7, 1, 6), cex.axis=0.7, pos=0)
arrows(0, 0.2, 2286389, 0.2, length=0, col="maroon", lwd=1.5)
points(c(1743190, 1747664), c(0.2, 0.2), type="p", cex=1.5, pch=19, col="maroon")

mtext("ChrR", side=3, adj=0.05, cex=0.85)
axis(1, cex.axis=0.8, at = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000), labels=c("0", "500", "1000", "1500", "2000", "2500", "3000"), cex.axis=0.8)
mtext("Chromosome position (bp x 1000)", side=1, outer=TRUE, line=1)
dev.off()
