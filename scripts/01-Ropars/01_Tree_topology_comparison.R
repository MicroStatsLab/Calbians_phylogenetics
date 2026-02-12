library(readxl)
library(TreeDist)
library("ape","PCPS")
library(phangorn)

unphased <- ape::read.tree("data_in/01-Ropars/Original_ropars_1_clade13.txt")
phased <- ape::read.tree("data_in/01-Ropars/Haplotypo_1_clade13.txt")
MLST_unphased <- ape::read.tree("data_in/01-Ropars/Ropars_MLST_1_clade13.txt") 

#resolve polytomies
MLST_unphased <- multi2di(MLST_unphased)

RF.dist(phased, unphased, check.labels=TRUE, rooted=FALSE,normalize=T)
#0.6551724
KF.dist(phased, unphased, check.labels=TRUE, rooted=FALSE)
#0.2801507

RF.dist(phased, MLST_unphased, check.labels=TRUE, rooted=FALSE,normalize=T)
#0.862069
KF.dist(phased, MLST_unphased, check.labels=TRUE, rooted=FALSE)
# 0.0125206

RF.dist(unphased, MLST_unphased, check.labels=TRUE, rooted=FALSE,normalize=T)
#0.862069
KF.dist(unphased, MLST_unphased, check.labels=TRUE, rooted=FALSE)
#0.2859665




TreeDistance(phased, MLST_unphased)

SharedPhylogeneticInfo(
  phased,
  unphased,
  normalize = T,
  reportMatching = T,
  diag = TRUE
)


SharedPhylogeneticInfo(
  phased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE,
  diag = TRUE
)

SharedPhylogeneticInfo(
  unphased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE,
  diag = TRUE
)

PhylogeneticInfoDistance(
  unphased,
  phased,
  normalize = T,
  reportMatching = FALSE
)

PhylogeneticInfoDistance(
  phased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)

PhylogeneticInfoDistance(
  unphased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)



DifferentPhylogeneticInfo(
  phased,
  unphased,
  normalize = F,
  reportMatching = FALSE
)

DifferentPhylogeneticInfo(
  phased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)



DifferentPhylogeneticInfo(
  unphased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)


ClusteringInfoDistance(
  phased,
  unphased,
  normalize = T,
  reportMatching = FALSE
)


ClusteringInfoDistance(
  phased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)


ClusteringInfoDistance(
  unphased,
  MLST_unphased,
  normalize = T,
  reportMatching = FALSE
)



DifferentPhylogeneticInfo(
  phased,
  unphased)
VisualizeMatching(DifferentPhylogeneticInfo,phased,
                  unphased,edge.cex =  0.4,show.tip.label = F,edge.width = 0.5)


MutualClusteringInfo(
  phased,
  unphased,
  normalize = T,
  reportMatching = FALSE,
  diag = TRUE
)

