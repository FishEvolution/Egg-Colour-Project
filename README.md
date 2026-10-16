# Egg-Colour-Project

Data Overview.xlsx - Contains all the original data for Egg Colour and environmental and parental care data.

Egg Colour Project.R - Overall R script file for all analyses. Separate scripts for each individual analyses are found within their respect folders.

Tree File - Due to size it is obtainable from Rabosky et al. 2018 (Rabosky DL, Chang J, Cowman PF, Sallan L, Friedman M, Kaschner K, Garilao C, Near TJ, Coll M, Alfaro ME. An inverse latitudinal gradient in speciation rate for marine fishes. Nature. 2018 Jul;559(7714):392-5.)

Below is a description of how each part of the methods were carried out including file names used in each analysis.

For all parts of the project the end of the file represents which part of the project it is used in:

No "1" - Data with the egg colours that match the colour of the adult were included.
1 - Data with the egg colours that match the colour of the adult were excluded.

A) Data Generation

1. Create a separate folder entitled maps.
2. Use the RangesMaps.csv file that contains the names of the .shp files needed to collect IUCN data.
3. Collect environmental data for each species in the RangeMaps.csv using the loop code.
4. Save the results in an Answers.csv file.
5. Repeat for all environmental variables (Temperature, Oxygen and Chlorophyll).
6. Repeat for AquaMaps using the RangeMapsCoords.csv file.

B) BayesTraits

1. Copy species names and egg colour pigment data (coded as N = Transparent (No pigment), C = Carotenoid (Red-Pink and Yellow-Amber-Orange), G = Green, M = Melanin (Black-Brown-Grey) and W = Guanine (White)) from FishEggs.csv / FishEggs1.csv for speices with egg colour data into a BayesTraits.csv and .txt file.
2. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape".
3. Produce a Pruned Tree of only species in the BayesTrait.csv file using the package "ape" and export it.
4. Remove all text from after the word #NEXUS until the BEGIN TREES; line from the .trees file.
5. Place the .txt file and .trees file in the BayesTraits folder.
6. Open the command window and the set directory to the BayesTraits folder using the cd command.
7. Run the following code: BayesTraitsV4.exe Tree.trees EggColour.txt
8. Choose the following options:
	1 (For the multistate model)
	2 (For MCMC)
	RevJump exp 10 (As there are many different parameters to test)
	EqualTrees 10000 (If one tree is being focused on)
	run

C) Ancestral State Reconstruction and Phylogenetic Signal

1. Create two ancestral state files with no headings, called AncestralState.csv which contains all the species names and egg colour data but with numerical categories (coded as 1 = Transparent, 2 = Carotenoid, 3 = Green, 4 = Guanine, 5 = Melanin).
2. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
3. Prune the Tree to only include species found in the GroupNames.csv file using the package "phangorn".
4. Create a maximum likelihood consensus tree using BEAST.
5. Ensure the tree is ultrametric.
6. Determine phylogenetic signal for traits using the package "caper" using the Fritz and Purvis D statistic.
7. Use FitDiscrete using the package "geiger" to test the different evolutionary models.
8. Perform an ancestral state reconstruction using the above models.

D) Mapping

1. Create a list of the five colours being mapped (Transparent, Carotenoid, Green, Guanine, Melanin).
2. For each colour, for all species, create a map of all points of occurrence, with each species shown in a different colour, using sp_occurrence in geodata and.
1. Create a ggplot beeswarm plot of environmental variables (y) against pigment (x).

E) Generalised Linear Mixed Models

1. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
2. Import FishEggs.csv and remove any rows that do not contain data for the variable being tested (Temperature, Oxygen Concentration, Chlorophyll Concentration, Latitude, Post Spawning Parental Care Method and Place of Development separately).
3. Run 5 separate analyses for each colour, and repeating this for data that includes and excludes eggs that match the colour of the parent. Note: Melanin has no eggs that do not match the colour of the adult and thus are not included in the additional analyses.
4. Perform a MCMCglmm analysis on the data and save the model as "Colour[NameOfColour][Variable].Rdata"

F) Path Analysis

1. Import FishEgg.csv.
2. Perform a PCA on the environmental variables (Temperature, Oxygen, Chlorophyll and Depth).
3. Write the results as a PCAScores.csv file and import these values into the FishEggs.csv file.
4. Load in the ConsensusTree.nex file.
5. For each pigment create 8 models. The first four: 1) all variables, 2) all variables and an interaction from environment to parental care and development, 3) all variables and an interaction from parental care to development and 4) all variables and all additional interactions from 2 and 3. The second four are the same but the predictor variables are restricted to those that were found to be associated with that egg pigment in our phylogenetic logistic regressions. Melanin, which had no predictor variables will only have four models.
6. Perform a path analysis using the phylo_path command.

G) Phylogenetic Signal

Import the PhyloSignal.csv file that contains all species and their associated pigment, and a consensus tree.
2. Create a column for each pigment with a 0 for pigment not present and a 1 for pigment present.
3. Perform a PhyloD analysis on each colour individually.
