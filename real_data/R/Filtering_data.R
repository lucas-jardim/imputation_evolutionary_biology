library("ape")
library("phytools")
library("taxadb")
library("dplyr")


# Download taxadb database ------------------------------------------------

td_create("gbif")


# Import phylogeny and phenotypic data ------------------------------------

phy <- ape::read.nexus("./real_data/data/Upham_phylogeny.nexus")
data <- read.csv("./real_data/data/primate_data.csv", header = TRUE)


# Clean names ---------------------------------------------------

phy_primates <- keep.tip(phy, grep("PRIMATES", phy$tip.label, value = TRUE))
phy_primates$tip.label <- gsub("^([A-Za-z]+_[a-z]+)_.*$", "\\1", phy_primates$tip.label) %>% gsub("_", " ", .)
data$Taxon <- gsub("_", " ", data$Taxon)   
data$Taxon[grep("Lagothrix lagotricha", data$Taxon)] <- "Lagothrix lagothricha"
phy_primates$tip.label[grep("Lagothrix lagotricha", phy_primates$tip.label)] <- "Lagothrix lagothricha"

# Cleaning duplicates -------------------------------------------

data$Taxon <- gsub("^([A-Za-z]+ [a-z]+) [a-z]+$", "\\1", data$Taxon)
data[grep("Cercocebus torquatus", data$Taxon)[1], 2:3] <- apply(data[grep("Cercocebus torquatus", data$Taxon), 2:3], 2, mean)
data <- data[-grep("Cercocebus torquatus", data$Taxon)[2], ]


# Checking for valid names -----------------------------------------------------
# Data --------------------------------------------------------------------

not_in_phylogeny <- data$Taxon[!(data$Taxon %in% phy_primates$tip.label)] 
data_db <- taxadb::filter_name(not_in_phylogeny, provider = "gbif")
accepted <- taxadb::filter_id(data_db$acceptedNameUsageID, "gbif")
accepted <- accepted[match(not_in_phylogeny, data_db[match(accepted$input, data_db$acceptedNameUsageID), "input"][[1]]), ] 
data[match(not_in_phylogeny, data$Taxon), "Taxon"] <- gsub("^([A-Za-z]+ [a-z]+) [a-z]+$", "\\1", accepted$scientificName) 

# Phylogeny ----------------------------------------------------------------

not_in_data <- phy_primates$tip.label[!(phy_primates$tip.label %in% data$Taxon)]
not_in_data <- not_in_data[!(not_in_data %in% c("Archaeolemur majori", "Homo denisova"))]
phylo_db <- taxadb::filter_name(not_in_data, "gbif")
phylo_db <- phylo_db[!is.na(phylo_db$scientificName), ]
phy_accepted <- filter_id(phylo_db$acceptedNameUsageID, "gbif")
phy_accepted <- phy_accepted[match(not_in_data, phylo_db[match(phy_accepted$input, phylo_db$acceptedNameUsageID), "input"][[1]]), ]  
phy_primates$tip.label[match(not_in_data, phy_primates$tip.label)] <-  gsub("^([A-Za-z]+ [a-z]+) [a-z]+$", "\\1", phy_accepted$scientificName) 


# Adding tip into phylogeny -----------------------------------------------

# Downloading reference trees ---------------------------------------------

otl_link <- c("https://api.opentreeoflife.org/v3/study/pg_2656.tre") #https://doi.org/10.1371/journal.pone.0049521
reference_trees <- read.newick(otl_link)[[2]]


# Imputing tips -----------------------------------------------------------

missing_taxa <- data$Taxon[!(data$Taxon %in% phy_primates$tip.label)]
phy_primates$tip.label[grep("Chiropotes utahickae", phy_primates$tip.label)] <- "Chiropotes satanas"
Cercocebus_agilis_branch <- phy_primates$edge.length[phy_primates$edge[, 2] == grep("Cercocebus agilis", phy_primates$tip.label)]
Cercocebus_galeritus_branch <- reference_trees$edge.length[reference_trees$edge[, 2]  == grep("Cercocebus galeritus", reference_trees$tip.label)]
phy_primates <- bind.tip(phy_primates, 
                          tip.label = "Cercocebus galeritus", 
                          where = grep("Cercocebus agilis", phy_primates$tip.label),
                          edge.length = Cercocebus_galeritus_branch,
                          position = Cercocebus_galeritus_branch)

Macaca_arctoides_branch <- phy_primates$edge.length[phy_primates$edge[, 2] == grep("Macaca arctoides", phy_primates$tip.label)]
Macaca_sinica_branch <- reference_trees$edge.length[reference_trees$edge[, 2]  == grep("Macaca sinica", reference_trees$tip.label)]
phy_primates <- bind.tip(phy_primates, 
                         tip.label = "Macaca sinica", 
                         where = grep("Macaca arctoides", phy_primates$tip.label),
                         edge.length = Macaca_sinica_branch,
                         position = Macaca_sinica_branch)

phy_primates <- keep.tip(phy_primates, data$Taxon)
data$Taxon <- gsub(" ", "_", data$Taxon)
write.csv(data, "./real_data/data/analysis_data.csv", row.names = FALSE)
write.tree(phy_primates, "./real_data/data/analysis_tree.tree")

