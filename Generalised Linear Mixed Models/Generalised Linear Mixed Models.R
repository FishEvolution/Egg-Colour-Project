#### Generalised Linear Mixed Models ####
#### Setup ####

library(ape)
library(regclass)
library(phangorn)
library(MCMCglmm)
library(coda)
library(ggplot2)
library(ggbeeswarm)

#### Temperature Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceTemp)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceTemp)),]}

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

FishEggs$zTemp<-scale(FishEggs$MeanSurfaceTemp)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentTemp.Rdata")
  
}

save(Final.Mod, file="ColourTransparentTemp.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidTemp.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidTemp.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninTemp.Rdata")
  
}

save(Final.Mod, file="ColourMelaninTemp.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineTemp.Rdata")
  
}

save(Final.Mod, file="ColourGuanineTemp.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenTemp.Rdata")
  
}

save(Final.Mod, file="ColourGreenTemp.Rdata")

summary(Final.Mod)

#### Temperature Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceTemp)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceTemp)),]}

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

FishEggs$zTemp<-scale(FishEggs$MeanSurfaceTemp)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentTemp1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentTemp1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidTemp1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidTemp1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineTemp1.Rdata")
  
}

save(Final.Mod, file="ColourGuanineTemp1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zTemp,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zTemp,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenTemp1.Rdata")
  
}

save(Final.Mod, file="ColourGreenTemp1.Rdata")

#### Depth Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$AverageDepth)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$AverageDepth)),]}

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

FishEggs$zDepth<-scale(FishEggs$AverageDepth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentDepth.Rdata")
  
}

save(Final.Mod, file="ColourTransparentDepth.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidDepth.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidDepth.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninDepth.Rdata")
  
}

save(Final.Mod, file="ColourMelaninDepth.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineDepth.Rdata")
  
}

save(Final.Mod, file="ColourGuanineDepth.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenDepth.Rdata")
  
}

save(Final.Mod, file="ColourGreenDepth.Rdata")

summary(Final.Mod)

#### Depth Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$AverageDepth)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$AverageDepth)),]}

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

FishEggs$zDepth<-scale(FishEggs$AverageDepth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentDepth1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentDepth1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentDepth1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidDepth1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidDepth1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidDepth1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineDepth1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineDepth1.Rdata")
  
}

save(Final.Mod, file="ColourGuanineDepth1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zDepth,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenDepth1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zDepth,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenDepth1.Rdata")
  
}

save(Final.Mod, file="ColourGreenDepth1.Rdata")

#### Oxygen Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceOxygen)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceOxygen)),]}

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

FishEggs$zOxygen<-scale(FishEggs$MeanSurfaceOxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentOxygen.Rdata")
  
}

save(Final.Mod, file="ColourTransparentOxygen.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidOxygen.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidOxygen.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninOxygen.Rdata")
  
}

save(Final.Mod, file="ColourMelaninOxygen.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineOxygen.Rdata")
  
}

save(Final.Mod, file="ColourGuanineOxygen.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenOxygen.Rdata")
  
}

save(Final.Mod, file="ColourGreenOxygen.Rdata")

summary(Final.Mod)

#### Oxygen Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceOxygen)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceOxygen)),]}

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

FishEggs$zOxygen<-scale(FishEggs$MeanSurfaceOxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentOxygen1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentOxygen1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidOxygen1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidOxygen1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineOxygen1.Rdata")
  
}

save(Final.Mod, file="ColourGuanineOxygen1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zOxygen,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zOxygen,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenOxygen1.Rdata")
  
}

save(Final.Mod, file="ColourGreenOxygen1.Rdata")
#### Chlorophyll Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceChlorophyll)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceChlorophyll)),]}

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

FishEggs$zChlorophyll<-scale(FishEggs$MeanSurfaceChlorophyll)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentChlorophyll.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentChlorophyll.Rdata")
  
}

save(Final.Mod, file="ColourTransparentChlorophyll.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidChlorophyll.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidChlorophyll.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidChlorophyll.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninChlorophyll.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninChlorophyll.Rdata")
  
}

save(Final.Mod, file="ColourMelaninChlorophyll.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineChlorophyll.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineChlorophyll.Rdata")
  
}

save(Final.Mod, file="ColourGuanineChlorophyll.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenChlorophyll.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenChlorophyll.Rdata")
  
}

save(Final.Mod, file="ColourGreenChlorophyll.Rdata")

summary(Final.Mod)

#### Chlorophyll Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$MeanSurfaceChlorophyll)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceChlorophyll)),]}

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

FishEggs$zChlorophyll<-scale(FishEggs$MeanSurfaceChlorophyll)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentChlorophyll1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentChlorophyll1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentChlorophyll1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidChlorophyll1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidChlorophyll1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidChlorophyll1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineChlorophyll1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineChlorophyll1.Rdata")
  
}

save(Final.Mod, file="ColourGuanineChlorophyll1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zChlorophyll,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zChlorophyll,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenChlorophyll1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zChlorophyll,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenChlorophyll1.Rdata")
  
}

save(Final.Mod, file="ColourGreenChlorophyll1.Rdata")
#### Positive Latitude Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$PositiveLatitude)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$PositiveLatitude)),]}

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

FishEggs$zPosLat<-scale(FishEggs$PositiveLatitude)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPosLat.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPosLat.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPosLat.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPosLat.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPosLat.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPosLat.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninPosLat.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninPosLat.Rdata")
  
}

save(Final.Mod, file="ColourMelaninPosLat.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePosLat.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePosLat.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePosLat.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPosLat.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPosLat.Rdata")
  
}

save(Final.Mod, file="ColourGreenPosLat.Rdata")

summary(Final.Mod)

#### Positive Latitude Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$PositiveLatitude)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$PositiveLatitude)),]}

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

FishEggs$zPosLat<-scale(FishEggs$PositiveLatitude)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPosLat1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPosLat1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPosLat1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPosLat1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPosLat1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPosLat1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePosLat1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePosLat1.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePosLat1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zPosLat,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPosLat,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPosLat1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPosLat,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPosLat1.Rdata")
  
}

save(Final.Mod, file="ColourGreenPosLat1.Rdata")


#### Parental Post Care ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$PostCare)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$PostCare)),]}

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

FishEggs$zPostCare<-as.factor(FishEggs$PostCare)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPostCare.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPostCare.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPostCare.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPostCare.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPostCare.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPostCare.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninPostCare.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninPostCare.Rdata")
  
}

save(Final.Mod, file="ColourMelaninPostCare.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePostCare.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePostCare.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePostCare.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPostCare.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPostCare.Rdata")
  
}

save(Final.Mod, file="ColourGreenPostCare.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

#### Parental Care No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$PostCare)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$PostCare)),]}

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

FishEggs$zPostCare<-as.factor(FishEggs$PostCare)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPostCare1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPostCare1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPostCare1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPostCare1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPostCare1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPostCare1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePostCare1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePostCare1.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePostCare1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPostCare,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPostCare,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPostCare1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPostCare,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPostCare1.Rdata")
  
}

save(Final.Mod, file="ColourGreenPostCare1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zPostCare1"] - Final.Mod$Sol[,"zPostCare2"]

HPDinterval(Difference)

#### Development Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$Development)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$Development)),]}

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

FishEggs$zDevelopment<-as.factor(FishEggs$DevelopmentScore)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentDevelopment.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentDevelopment.Rdata")
  
}

save(Final.Mod, file="ColourTransparentDevelopment.Rdata")

summary(Final.Mod)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidDevelopment.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidDevelopment.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidDevelopment.Rdata")

summary(Final.Mod)

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninDevelopment.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninDevelopment.Rdata")
  
}

save(Final.Mod, file="ColourMelaninDevelopment.Rdata")

summary(Final.Mod)

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineDevelopment.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineDevelopment.Rdata")
  
}

save(Final.Mod, file="ColourGuanineDevelopment.Rdata")

summary(Final.Mod)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenDevelopment.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenDevelopment.Rdata")
  
}

save(Final.Mod, file="ColourGreenDevelopment.Rdata")

summary(Final.Mod)
#### Development Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]
if(sum(is.na(FishEggs$Development)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$Development)),]}

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

FishEggs$zDevelopment<-as.factor(FishEggs$DevelopmentScore)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentDevelopment1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentDevelopment1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentDevelopment1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zDevelopment1"] - Final.Mod$Sol[,"zDevelopment2"]

HPDinterval(Difference)

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidDevelopment1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidDevelopment1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidDevelopment1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zDevelopment1"] - Final.Mod$Sol[,"zDevelopment2"]

HPDinterval(Difference)
# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuanineDevelopment1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuanineDevelopment1.Rdata")
  
}

save(Final.Mod, file="ColourGuanineDevelopment1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zDevelopment1"] - Final.Mod$Sol[,"zDevelopment2"]

HPDinterval(Difference)

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zDevelopment,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zDevelopment,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenDevelopment1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zDevelopment,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenDevelopment1.Rdata")
  
}

save(Final.Mod, file="ColourGreenDevelopment1.Rdata")

summary(Final.Mod)

Difference <- Final.Mod$Sol[,"zDevelopment1"] - Final.Mod$Sol[,"zDevelopment2"]

HPDinterval(Difference)

#### Temperature + Oxygen PC Model ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

PCA <- FishEggs[c(9:10)]

if(sum(is.na(PCA$MeanSurfaceTemp)>0)){PCA<-PCA[-which(is.na(PCA$MeanSurfaceTemp)),]}

Data.PCA <- princomp(PCA, cor = TRUE)

summary(Data.PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]

if(sum(is.na(FishEggs$MeanSurfaceTemp)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceTemp)),]}

FishEggs$zPC1<-scale(FishEggs$PC1)
FishEggs$zPC2<-scale(FishEggs$PC2)

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

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPCtempOxy.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPCtempOxy.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPCtempOxy.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPCtempOxy.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPCtempOxy.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPCtempOxy.Rdata")

# Melanin

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Melanin)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourMelaninPCtempOxy.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Melanin)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourMelaninPCtempOxy.Rdata")
  
}

save(Final.Mod, file="ColourMelaninPCtempOxy.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePCtempOxy.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePCtempOxy.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePCtempOxy.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPCtempOxy.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPCtempOxy.Rdata")
  
}

save(Final.Mod, file="ColourGreenPCtempOxy.Rdata")

#### Temperature + Oxygen PC Model No Matching ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

FishEggs <- read.csv("FishEggs1.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

FishEggs <- FishEggs[FishEggs$Pigment != "N/A", ]

if(sum(is.na(FishEggs$MeanSurfaceTemp)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceTemp)),]}

FishEggs$zPC1<-scale(FishEggs$PC1)
FishEggs$zPC2<-scale(FishEggs$PC2)

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

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, FishEggs$Species))
Tree <- T100[[1]]

FishEggs$Species %in% Tree$tip.label
FishEggs<-FishEggs[(FishEggs$Species %in% Tree$tip.label),]

# Transparent

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Transparent)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourTransparentPCTempOxy1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Transparent)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourTransparentPCTempOxy1.Rdata")
  
}

save(Final.Mod, file="ColourTransparentPCTempOxy1.Rdata")

# Carotenoid

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2, 
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Carotenoid)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourCarotenoidPCTempOxy1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Carotenoid)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourCarotenoidPCTempOxy1.Rdata")
  
}

save(Final.Mod, file="ColourCarotenoidPCTempOxy1.Rdata")

# Guanine

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Guanine)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGuaninePCTempOxy1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Guanine)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGuaninePCTempOxy1.Rdata")
  
}

save(Final.Mod, file="ColourGuaninePCTempOxy1.Rdata")

# Green

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = FishEggs, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Mod<-MCMCglmm(as.factor(Green)~zPC1 + zPC2,
              random = ~Species, 
              ginverse = list(Species = animalA), 
              prior = gelmanprior, 
              verbose = TRUE, 
              family = "categorical", 
              data = FishEggs,
              nitt = 11000,
              thin = 10,
              burnin = 1000,
              pl = TRUE, 
              pr = TRUE, 
              slice = TRUE) 

Final.Mod<-Mod
Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 

nsamp.l <- nrow(Mod$VCV)
start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l, "Species"]))

save(Final.Mod, file = "ColourGreenPCTempOxy1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Mod <- MCMCglmm(as.factor(Green)~zPC1 + zPC2,
                  random = ~Species, 
                  ginverse=list(Species = animalA), 
                  prior = gelmanprior, 
                  verbose = FALSE, 
                  family = "categorical",  
                  start = start1.l,
                  data = FishEggs,
                  nitt = 22000,
                  thin = 2000,
                  burnin = 2000,
                  pl = TRUE,
                  pr = TRUE,
                  slice = TRUE)
  
  print(i)
  
  Final.Mod$VCV[((i-1)*10+1):(i*10), ] <- Mod$VCV[1:10,] 
  Final.Mod$Sol[((i-1)*10+1):(i*10), ] <- Mod$Sol[1:10,] 
  Final.Mod$Liab[((i-1)*10+1):(i*10), ] <- Mod$Liab[1:10,] 
  
  nsamp.l <- nrow(Mod$VCV)
  start1.l = list(R = Mod$VCV[nsamp.l, "units"], G = list(G1 = Mod$VCV[nsamp.l,"Species"]))
  
  save(Final.Mod, file = "ColourGreenPCTempOxy1.Rdata")
  
}

save(Final.Mod, file="ColourGreenPCTempOxy1.Rdata")
#### Graph ####

FishEggs <- read.csv("FishEggs.csv")
attach(FishEggs)

if(sum(is.na(FishEggs$Pigment)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$Pigment)),]} 
if(sum(is.na(FishEggs$MeanSurfaceTemp)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$MeanSurfaceTemp)),]} 

ggplot(data = FishEggs) +
  aes(y = MeanSurfaceTemp, x = Pigment) +
  geom_point(size=5) +
  geom_beeswarm(cex = 0.1) +
  theme(text = element_text(size = 35),
        axis.text.x = element_text(hjust = 1)) +
  labs(y = "Mean Surface Temperature (\u00B0C)", x = "Pigment") +
  coord_flip()

by(FishEggs$MeanSurfaceTemp, FishEggs$Pigment, median)
by(FishEggs$MeanSurfaceTemp, FishEggs$Pigment, quantile, probs=c(0.25,0.75), na.rm=TRUE)
