#### Ancestral State Reconstruction ####
#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Ancestral State Reconstruction")

library(ape)
library(phangorn)
library(picante)
library(geiger)
library(geomorph)

#### Consensus Tree ####

Phylogeny <- read.tree("Actinopterygii.trees")

GroupNames <- read.csv("GroupNames.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Phylogeny[[1]]$tip.label)

for(k in 1:100)
{
  Phylogeny[[k]] <- keep.tip(Phylogeny[[k]], Phylogeny[[k]]$tip.label[match(GroupNames$V1, Phylogeny[[k]]$tip.label)])
  print(k)
}

write.nexus(Phylogeny, file="Phylogeny.nex")

ConsensusTree <- read.nexus("ConsensusTree.nex")

ConsensusTree <- nnls.tree(method = "ultrametric", cophenetic(ConsensusTree), ConsensusTree)

min(ConsensusTree$edge.length)

ConsensusTree$edge.length[ConsensusTree$edge.length<0.0001] <- 0.0001

plot(ConsensusTree)

#### Ancestral States ####

AncestralStates <- read.csv("AncestralStates.csv", row.names=1, header=FALSE)

fitDiscrete(ConsensusTree, AncestralStates, model = c("ER"))
fitDiscrete(ConsensusTree, AncestralStates, model = c("SYM"))
fitDiscrete(ConsensusTree, AncestralStates, model = c("ARD"))
fitDiscrete(ConsensusTree, AncestralStates, model = c("ARD"), transform = "delta")

# ARD No Delta

Fit <- ace(AncestralStates$V2, ConsensusTree, model = "ARD", type = "discrete")
Fit
round(Fit$lik.anc, 3)

palette(c("white", "orange", "green4", "grey", "black"))
palette()

plot.phylo(ConsensusTree, type = "fan", cex = 0.4, label.offset = 20, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStates$V2))],
               sort(unique(AncestralStates$V2)))
par(lwd = 0.1)
nodelabels(node = 1:ConsensusTree$Nnode+Ntip(ConsensusTree), pie = Fit$lik.anc, 
           piecol = cols, cex = 0.15)
tiplabels(pie=to.matrix(AncestralStates$V2, sort(unique(AncestralStates$V2))), 
          piecol = cols, cex = 0.15)