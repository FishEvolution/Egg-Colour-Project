#### Setup ####

library(ape)
library(phytools)
library(phangorn)
library(geiger)
library(RRphylo)
library(caper)
library(terra)
library(geodata)
library(MCMCglmm)
library(car)
library(ggplot2)
library(ggbeeswarm)
library(raster)
library(sf)
library(lavaan)
library(semPlot)
library(graph)
library(phylopath)

#### Data Collection ####

# Egg Colour URL Collection

EggColourCollection <- read.csv("EggColourCollection.csv")

N <- dim(EggColourCollection)[1]

answers <- matrix(NA, nrow=N, ncol=2)

for(i in 1:N){
  
  answers[i,1] <- i
  
  GenusName <- EggColourCollection$Genus[i]
  SpeciesName <- EggColourCollection$Species[i]
  Code1 <- EggColourCollection$SpecCode[i]
  Code2 <- EggColourCollection$StockCode[i]
  URL <- paste0("https://www.fishbase.se/Reproduction/FishEggInfoSummary.php?ID=",Code1,"&GenusName=",GenusName,"&SpeciesName=",SpeciesName,"&fc=1&StockCode=",Code2)
  
  answers[i,2] <- URL
  
  print(i)
}

write.csv(answers,file="URLS.csv")

urls <- read.csv("URLS.csv", header =  FALSE)

urls <- urls$V1

valid_url <- function(url_in,t=2){
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=T)[1])
  suppressWarnings(try(close.connection(con),silent=T))
  ifelse(is.null(check),TRUE,FALSE)
}

URL <- sapply(urls,valid_url)

write.csv(URL,file="URLExist.csv")

# Egg Colour Data

EggColourCollection <- read.csv("EggColourCollection.csv")

N <- dim(EggColourCollection)[1]

answers <- matrix(NA, nrow=N, ncol=2)

for(i in 1712:N){
  
  answers[i,1] <- i
  
  GenusName <- EggColourCollection$Genus[i]
  SpeciesName <- EggColourCollection$Species[i]
  Code1 <- EggColourCollection$SpecCode[i]
  Code2 <- EggColourCollection$StockCode[i]
  URL <- paste0("https://www.fishbase.se/Reproduction/FishEggInfoSummary.php?ID=",Code1,"&GenusName=",GenusName,"&SpeciesName=",SpeciesName,"&fc=1&StockCode=",Code2)
  
  Text <- read_html(URL)
  
  Data <- Text %>% html_elements("td.threefourth")
  
  Data <- Data[1]
  
  EggColourData <- data.frame(hrefs=as(Data, "character"))
  
  if(is_empty(EggColourData$hrefs) == 'FALSE') {
    answers[i,2] <- EggColourData$hrefs
  } else {
    answers[i, 2] <- NA_real_
  }
  
  print(i)
  Sys.sleep(0.5)
}

write.csv(answers,file="EggColourData.csv")

# Environmental Data

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/R/Maps")

# IUCN Files

x <- list.files()
x
Temperature <- x[225]
Oxygen <- x[181]
Chlorophyll <- x[41]

Temp <- raster(Temperature)
Oxy <- raster(Oxygen)
Chlor <- raster(Chlorophyll)

ListMaps<-read.csv("RangeMaps.csv", header=FALSE)
ListMaps<-as.vector(ListMaps)

N <- length(ListMaps)
TempAnswers <- rep(0,N)
OxyAnswers <- rep(0,N)
ChlorAnswers <- rep(0,N)

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Temp,Shape),Shape)
  answer<-mean(getValues(y2),na.rm=TRUE) # Take SD, Simulate from Normal Distribution. MCMCglmm - 10 values per species
  
  TempAnswers[i]<-mean(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(TempAnswers, file=paste("TempAnswers.csv"))

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Oxy,Shape),Shape)
  answer<-mean(getValues(y2),na.rm=TRUE)
  
  OxyAnswers[i]<-mean(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(OxyAnswers, file=paste("OxyAnswers.csv"))

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Chlor,Shape),Shape)
  answer<-mean(getValues(y2),na.rm=TRUE)
  
  ChlorAnswers[i]<-mean(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(ChlorAnswers, file=paste("ChlorAnswers.csv"))

# AquaMaps Files

x <- list.files()
x
Temperature <- x[213]
Oxygen <- x[166]
Chlorophyll <- x[31]

Temp <- raster(Temperature)
Oxy <- raster(Oxygen)
Chlor <- raster(Chlorophyll)

ListMapsCoords <- read.csv("RangeMapsCoords.csv", header=FALSE)
ListMapsCoords <- c(ListMapsCoords$V1)

N <- length(ListMapsCoords)

TempAnswersCoords <- rep(0,N)
OxyAnswersCoords <- rep(0,N)
ChlorAnswersCoords <- rep(0,N)

for(i in 1:N){
  Species <- read.csv(ListMapsCoords[i])
  Species <- Species[Species$Probability >= 0.5, ]
  
  Coordinates <- data.frame(lon=Species$Center.Long,
                            lat = Species$Center.Lat)
  coordinates(Coordinates)<-c("lon","lat")
  
  Values <- extract(x = Temp, y = Coordinates)
  Values <- Values[!is.na(Values)]
  
  answer <- mean(Values)
  
  TempAnswersCoords[i] <- mean(Values)
  print(i)
}

write.csv(TempAnswersCoords, file=paste("TempAnswersCoords.csv"))

for(i in 1:N){
  Species <- read.csv(ListMapsCoords[i])
  Species <- Species[Species$Probability >= 0.5, ]
  
  Coordinates <- data.frame(lon=Species$Center.Long,
                            lat = Species$Center.Lat)
  coordinates(Coordinates)<-c("lon","lat")
  
  Values <- extract(x = Oxy, y = Coordinates)
  Values <- Values[!is.na(Values)]
  
  OxyAnswersCoords[i] <- mean(Values)
  print(i)
}

write.csv(OxyAnswersCoords, file=paste("OxyAnswersCoords.csv"))

for(i in 1:N){
  Species <- read.csv(ListMapsCoords[i])
  Species <- Species[Species$Probability >= 0.5, ]
  
  Coordinates <- data.frame(lon=Species$Center.Long,
                            lat = Species$Center.Lat)
  coordinates(Coordinates)<-c("lon","lat")
  
  Values <- extract(x = Chlor, y = Coordinates)
  Values <- Values[!is.na(Values)]
  
  ChlorAnswersCoords[i] <- mean(Values)
  print(i)
}

write.csv(ChlorAnswersCoords, file=paste("ChlorAnswersCoords.csv"))

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

#### Ancestral State Reconstruction ####

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
#### Phylogenetic Signal ####

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

#### Generate Maps ####

Dir <- "C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/Fish_Raster_Stack.tif"
Map <- rast(Dir) 

MapNames <- read.csv("MapNames.csv")
names(Map) <- gsub(" ","_",MapNames$Name)

FishData <- read.csv("FishEggColourForMapping.csv")

FishMap <- FishData[match(names(Map), FishData$Species), ] 

# None

Fdisp_fun <- function(cell_values) {
  Present <- which(cell_values == 1)
  if (length(Present) == 0) return(0)
  spp <- FishMap$Species[Present]
  spp <- spp[spp %in% FishData$Species]
  
  Traits <- FishData[which(FishData$Species %in% spp),c("None")]
  
  Prop <- sum(Traits)/length(spp)
  
  return(Prop)
}

Fdisp_raster <- terra::app(Map, Fdisp_fun)

saveRDS(Fdisp_raster, file = "None.rds")

Fdisp_raster <- readRDS("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/None.rds")

None <- Fdisp_raster
None[is.na(Fdisp_raster)] <- 0

None[None == 0] <- NA

world <- ne_countries(scale = "medium", returnclass = "sf")
world_eq <- st_transform(world, crs = 6933)
plot(None, col = map.pal("plasma", 100), legend = TRUE, axes = TRUE)
plot(vect(world_eq), add = TRUE, border = "black", col = NA)

# Carotenoid

Fdisp_fun <- function(cell_values) {
  Present <- which(cell_values == 1)
  if (length(Present) == 0) return(0)
  spp <- FishMap$Species[Present]
  spp <- spp[spp %in% FishData$Species]
  
  Traits <- FishData[which(FishData$Species %in% spp),c("Carotenoid")]
  
  Prop <- sum(Traits)/length(spp)
  
  return(Prop)
}

Fdisp_raster <- terra::app(Map, Fdisp_fun)

saveRDS(Fdisp_raster, file = "Carotenoid.rds")

Fdisp_raster <- readRDS("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/Carotenoid.rds")

Carotenoid <- Fdisp_raster
Carotenoid[is.na(Fdisp_raster)] <- 0

Carotenoid[Carotenoid == 0] <- NA

world <- ne_countries(scale = "medium", returnclass = "sf")
world_eq <- st_transform(world, crs = 6933)
plot(Carotenoid, col = map.pal("plasma", 100), legend = TRUE, axes = TRUE)
plot(vect(world_eq), add = TRUE, border = "black", col = NA)

# Guanine

Fdisp_fun <- function(cell_values) {
  Present <- which(cell_values == 1)
  if (length(Present) == 0) return(0)
  spp <- FishMap$Species[Present]
  spp <- spp[spp %in% FishData$Species]
  
  Traits <- FishData[which(FishData$Species %in% spp),c("Guanine")]
  
  Prop <- sum(Traits)/length(spp)
  
  return(Prop)
}

Fdisp_raster <- terra::app(Map, Fdisp_fun)

saveRDS(Fdisp_raster, file = "Guanine.rds")

Fdisp_raster <- readRDS("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/Guanine.rds")

Guanine <- Fdisp_raster
Guanine[is.na(Fdisp_raster)] <- 0

Guanine[Guanine == 0] <- NA

world <- ne_countries(scale = "medium", returnclass = "sf")
world_eq <- st_transform(world, crs = 6933)
plot(Guanine, col = map.pal("plasma", 100), legend = TRUE, axes = TRUE)
plot(vect(world_eq), add = TRUE, border = "black", col = NA)

# Melanin

Fdisp_fun <- function(cell_values) {
  Present <- which(cell_values == 1)
  if (length(Present) == 0) return(0)
  spp <- FishMap$Species[Present]
  spp <- spp[spp %in% FishData$Species]
  
  Traits <- FishData[which(FishData$Species %in% spp),c("Melanin")]
  
  Prop <- sum(Traits)/length(spp)
  
  return(Prop)
}

Fdisp_raster <- terra::app(Map, Fdisp_fun)

saveRDS(Fdisp_raster, file = "Melanin.rds")

Fdisp_raster <- readRDS("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/Melanin.rds")

Melanin <- Fdisp_raster
Melanin[is.na(Fdisp_raster)] <- 0

Melanin[Melanin == 0] <- NA

world <- ne_countries(scale = "medium", returnclass = "sf")
world_eq <- st_transform(world, crs = 6933)
plot(Melanin, col = map.pal("plasma", 100), legend = TRUE, axes = TRUE)
plot(vect(world_eq), add = TRUE, border = "black", col = NA)

# Green

Fdisp_fun <- function(cell_values) {
  Present <- which(cell_values == 1)
  if (length(Present) == 0) return(0)
  spp <- FishMap$Species[Present]
  spp <- spp[spp %in% FishData$Species]
  
  Traits <- FishData[which(FishData$Species %in% spp),c("Green")]
  
  Prop <- sum(Traits)/length(spp)
  
  return(Prop)
}

Fdisp_raster <- terra::app(Map, Fdisp_fun)

saveRDS(Fdisp_raster, file = "Green.rds")

Fdisp_raster <- readRDS("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping/Green.rds")

Green <- Fdisp_raster
Green[is.na(Fdisp_raster)] <- 0

Green[Green == 0] <- NA

world <- ne_countries(scale = "medium", returnclass = "sf")
world_eq <- st_transform(world, crs = 6933)
plot(Green, col = map.pal("plasma", 100), legend = TRUE, axes = TRUE)
plot(vect(world_eq), add = TRUE, border = "black", col = NA)

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
#### Path Analysis ####

FishEggs <- read.csv("FishEggs.csv")
colnames(FishEggs)[3] <- "Taxon-Family" 

PCA <- FishEggs[c(7:10)]

Data.PCA <- princomp(PCA, cor = TRUE)

summary(Data.PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

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
#### Graphs ####

FishEggs <- read.csv("FishEggs.csv")
attach(FishEggs)

if(sum(is.na(FishEggs$EggColourCategories)>0)){FishEggs<-FishEggs[-which(is.na(FishEggs$EggColourCategories)),]} 

ggplot(data = FishEggs) +
  aes(y = AverageDepth, x = Pigment) +
  geom_point(size=5) +
  geom_beeswarm(cex = 0.1) +
  theme(text = element_text(size = 35),
        axis.text.x = element_text(hjust = 1)) +
  labs(y = "Mean Surface Temperature (\u00B0C)", x = "Pigment") +
  coord_flip()

