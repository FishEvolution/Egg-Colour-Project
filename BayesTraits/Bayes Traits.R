#### Bayes Traits ####
#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/BayesTraits")

library(ape)

Phylogeny <- read.tree("Actinopterygii.trees")

#### BayesTraits Matching ####

Tree <- Phylogeny

BayesTraits <- read.csv("BayesTraits.csv", header=FALSE)

GroupNames <- read.csv("GroupNames.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Tree[[1]]$tip.label)

for(k in 1:100)
{
  Tree[[k]] <- keep.tip(Tree[[k]], Tree[[k]]$tip.label[match(GroupNames$V1, Tree[[k]]$tip.label)])
  print(k)
}

write.nexus(Tree, file="Tree.trees")

#### BayesTraits No Matching ####

Tree1 <- Phylogeny

BayesTraits1 <- read.csv("BayesTraits1.csv", header=FALSE)

GroupNames1 <- read.csv("GroupNames1.csv", header = F)

GroupNames1 <- subset(GroupNames1, V1 %in% Tree1[[1]]$tip.label)

for(k in 1:100)
{
  Tree1[[k]] <- keep.tip(Tree1[[k]], Tree1[[k]]$tip.label[match(BayesTraits1$V1, Tree1[[k]]$tip.label)])
  print(k)
}

write.nexus(Tree1, file = "Tree1.trees")
