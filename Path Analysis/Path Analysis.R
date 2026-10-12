#### Path Anaylsis ####
#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Path Analysis")

library(ape)
library(phangorn)
library(phylopath)

#### Consensus Trees ####

Phylogeny <- read.tree("Actinopterygii.trees")

GroupNames <- read.csv("GroupNames.csv", header = TRUE)

GroupNames <- subset(GroupNames, Species %in% Phylogeny[[1]]$tip.label)

for(k in 1:100)
{
  Phylogeny[[k]] <- keep.tip(Phylogeny[[k]], Phylogeny[[k]]$tip.label[match(GroupNames$Species, Phylogeny[[k]]$tip.label)])
  print(k)
}

write.nexus(Phylogeny, file="Phylogeny.nex")

ConsensusTree <- read.nexus("ConsensusTree.nex")

ConsensusTree <- nnls.tree(method = "ultrametric", cophenetic(ConsensusTree), ConsensusTree)

min(ConsensusTree$edge.length)

ConsensusTree$edge.length[ConsensusTree$edge.length<0.0001] <- 0.0001

plot(ConsensusTree)

#### PCA ####

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

PCA <- FishEggs[c(7:10)]

Data.PCA <- princomp(PCA, cor = TRUE)

summary(Data.PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

#### Analyses ####

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 
rownames(FishEggs) <- FishEggs$Species

FishEggs <- subset(FishEggs, FishEggs$Species %in% ConsensusTree$tip.label)

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

FishEggs$zPC1<-scale(FishEggs$PC1)
FishEggs$zPC2<-scale(FishEggs$PC2)
FishEggs$zPostCare<-as.factor(FishEggs$PostCare)
FishEggs$zDevelopment<-as.factor(FishEggs$Development)

# Transparent

ModelsT <- define_model_set(
  one = c(Transparent ~ zPC1 + zPC2 + zDevelopment + zPostCare),
  two = c(Transparent ~ zPC1 + zPC2 + zDevelopment + zPostCare,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  three = c(Transparent ~ zPC1 + zPC2 + zDevelopment + zPostCare,
            zDevelopment ~ zPostCare),
  four = c(Transparent ~ zPC1 + zPC2 + zDevelopment + zPostCare,
           zDevelopment ~ zPostCare,
           zDevelopment ~ zPC1 + zPC2,
           zPostCare ~ zPC1 + zPC2),
  five = c(Transparent ~ zPostCare + zDevelopment),
  six = c(Transparent ~ zPostCare + zDevelopment,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  seven = c(Transparent ~ zPostCare + zDevelopment,
            zDevelopment ~ zPostCare),
  eight = c(Transparent ~ zPostCare + zDevelopment,
            zDevelopment ~ zPostCare,
            zDevelopment ~ zPC1 + zPC2,
            zPostCare ~ zPC1 + zPC2))

plot_model_set(ModelsT)

resultT <- phylo_path(ModelsT, data = FishEggs, tree = ConsensusTree, model = 'lambda', method = "logistic_MPLE", btol = 40)

(T <- summary(resultT))

plot(T, main = "Transparent", cex.axis=100)

(best_modelT <- best(resultT))

plot(best_modelT, cex=100)

# Carotenoid

ModelsC <- define_model_set(
  one = c(Carotenoid ~ zPC1 + zPC2 + zDevelopment + zPostCare),
  two = c(Carotenoid ~ zPC1 + zPC2 + zDevelopment + zPostCare,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  three = c(Carotenoid ~ zPC1 + zPC2 + zDevelopment + zPostCare,
            zDevelopment ~ zPostCare),
  four = c(Carotenoid ~ zPC1 + zPC2 + zDevelopment + zPostCare,
           zDevelopment ~ zPostCare,
           zDevelopment ~ zPC1 + zPC2,
           zPostCare ~ zPC1 + zPC2),
  five = c(Carotenoid ~ zPostCare + zDevelopment),
  six = c(Carotenoid ~ zPostCare + zDevelopment,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  seven = c(Carotenoid ~ zPostCare + zDevelopment,
            zDevelopment ~ zPostCare),
  eight = c(Carotenoid ~ zPostCare + zDevelopment,
            zDevelopment ~ zPostCare,
            zDevelopment ~ zPC1 + zPC2,
            zPostCare ~ zPC1 + zPC2)
)

plot_model_set(ModelsC)

resultC <- phylo_path(ModelsC, data = FishEggs, tree = ConsensusTree, model = 'lambda', method = "logistic_MPLE", btol = 40)

(C <- summary(resultC))

plot(C)

(best_modelC <- best(resultC))

plot(best_modelC)

# Guanine

ModelsW <- define_model_set(
  one = c(Guanine ~ zPC1 + zPC2 + zDevelopment + zPostCare),
  two = c(Guanine ~ zPC1 + zPC2 + zDevelopment + zPostCare,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  three = c(Guanine ~ zPC1 + zPC2 + zDevelopment + zPostCare,
            zDevelopment ~ zPostCare),
  four = c(Guanine ~ zPC1 + zPC2 + zDevelopment + zPostCare,
           zDevelopment ~ zPostCare,
           zDevelopment ~ zPC1 + zPC2,
           zPostCare ~ zPC1 + zPC2),
  five = c(Guanine ~ zPC1 + zPC2 + zDevelopment),
  six = c(Guanine ~ zPC1 + zPC2 + zDevelopment,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  seven = c(Guanine ~ zPC1 + zPC2 + zDevelopment,
            zDevelopment ~ zPostCare),
  eight = c(Guanine ~ zPC1 + zPC2 + zDevelopment,
            zDevelopment ~ zPostCare,
            zDevelopment ~ zPC1 + zPC2,
            zPostCare ~ zPC1 + zPC2)
)

plot_model_set(ModelsW)

resultW <- phylo_path(ModelsW, data = FishEggs,tree = ConsensusTree, model = 'lambda', method = "logistic_MPLE", btol = 40)

(W <- summary(resultW))

plot(W)

(best_modelW <- best(resultW))

plot(best_modelW)

# Melanin - Just run the first four models

ModelsM <- define_model_set(
  one = c(Melanin ~ zPC1 + zPC2 + zDevelopment + zPostCare),
  two = c(Melanin ~ zPC1 + zPC2 + zDevelopment + zPostCare,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  three = c(Melanin ~ zPC1 + zPC2 + zDevelopment + zPostCare,
            zDevelopment ~ zPostCare),
  four = c(Melanin ~ zPC1 + zPC2 + zDevelopment + zPostCare,
           zDevelopment ~ zPostCare,
           zDevelopment ~ zPC1 + zPC2,
           zPostCare ~ zPC1 + zPC2)
)

plot_model_set(ModelsM)

resultM <- phylo_path(ModelsM, data = FishEggs, tree = ConsensusTree, model = 'lambda')

(M <- summary(resultM))

plot(M)

(best_modelM <- best(resultM))

plot(best_modelM)

# Green

ModelsG <- define_model_set(
  one = c(Green ~ zPC1 + zPC2 + zDevelopment + zPostCare),
  two = c(Green ~ zPC1 + zPC2 + zDevelopment + zPostCare,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  three = c(Green ~ zPC1 + zPC2 + zDevelopment + zPostCare,
            zDevelopment ~ zPostCare),
  four = c(Green ~ zPC1 + zPC2 + zDevelopment + zPostCare,
           zDevelopment ~ zPostCare,
           zDevelopment ~ zPC1 + zPC2,
           zPostCare ~ zPC1 + zPC2),
  five = c(Green ~ zDevelopment),
  six = c(Green ~ zDevelopment,
          zDevelopment ~ zPC1 + zPC2,
          zPostCare ~ zPC1 + zPC2),
  seven = c(Green ~ zDevelopment,
            zDevelopment ~ zPostCare),
  eight = c(Green ~ zDevelopment,
            zDevelopment ~ zPostCare,
            zDevelopment ~ zPC1 + zPC2,
            zPostCare ~ zPC1 + zPC2)
)

plot_model_set(ModelsG)

resultG <- phylo_path(ModelsG, data = FishEggs, tree = ConsensusTree, model = 'lambda')

(G <- summary(resultG))

plot(G)

(best_modelG <- best(resultG))

plot(best_modelG)