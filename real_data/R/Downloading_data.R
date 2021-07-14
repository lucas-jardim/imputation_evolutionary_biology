library("readxl")
library("dplyr")


# Download DeCasien data --------------------------------------------------

dir.create("real_data/data", recursive = TRUE)
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-017-0112/MediaObjects/41559_2017_BFs415590170112_MOESM250_ESM.xls",
              destfile = "./real_data/data/DeCasien_data.xls")


# Filtering brain and body mass data --------------------------------------

brain <- read_xls(path = "./real_data/data/DeCasien_data.xls", 
                  sheet = 3
                  ) %>% select(Taxon, `Final Brain Weight (g)`)
                               
                               
                               
body <- read_xls(path = "./real_data/data/DeCasien_data.xls", 
                 sheet = 4
                 )%>% select(Taxon, `Final Body Weight (g)`)


primate_data <- full_join(brain, body, by = "Taxon")

write.csv(primate_data, "./real_data/data/primate_data.csv", row.names = FALSE)


# Download phylogeny ------------------------------------------------------

download.file("https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre", 
              "./real_data/data/Upham_phylogeny.nexus")

