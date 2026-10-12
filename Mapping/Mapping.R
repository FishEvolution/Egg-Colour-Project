#### Mapping ####
#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Actinopteygii Egg Morphology/Data Repository/Mapping")

library(raster)
library(terra)
library(rnaturalearth)
library(ape)
library(fundiversity)

#### Generate Maps 2 ####

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