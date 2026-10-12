#### Phylogenetic Signal ####
#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Phylogenetic Signal")

library(ape)
library(caper)

#### PhyloD ####

ConsensusTree <- read.nexus("ConsensusTree.nex")

FishEggs <- read.csv("PhyloSignal.csv")

FishEggs$Transparent<-rep(0,dim(FishEggs)[1])
FishEggs$Transparent[which(FishEggs$Pigment=="None")]<-1

FishEggs$Carotenoid<-rep(0,dim(FishEggs)[1])
FishEggs$Carotenoid[which(FishEggs$Pigment=="Carotenoid")]<-1

FishEggs$Melanin<-rep(0,dim(FishEggs)[1])
FishEggs$Melanin[which(FishEggs$Pigment=="Melanin")]<-1

FishEggs$Guanine<-rep(0,dim(FishEggs)[1])
FishEggs$Guanine[which(FishEggs$Pigment=="Guanine")]<-1

FishEggs$Green<-rep(0,dim(FishEggs)[1])
FishEggs$Green[which(FishEggs$Pigment=="Green")]<-1

PhyloD <- comparative.data(ConsensusTree, FishEggs, Species)

TransparentPhyloD <- phylo.d(PhyloD, binvar=Transparent)
TransparentPhyloD

CarotenoidPhyloD <- phylo.d(PhyloD, binvar=Carotenoid)
CarotenoidPhyloD

MelaninPhyloD <- phylo.d(PhyloD, binvar=Melanin)
MelaninPhyloD

GuaninePhyloD <- phylo.d(PhyloD, binvar=Guanine)
GuaninePhyloD

GreenPhyloD <- phylo.d(PhyloD, binvar=Green)
GreenPhyloD
