# Loading packages --------------------------------------------------------

pacotes <- c("mice","geiger","ape","phytools","testthat", "phylolm","MPSEM","doParallel","letsR","vegan")
invisible(sapply(pacotes,
                 FUN= function(x) {if(!require(x,character.only=T)){install.packages(x,dep=T)
                   require(x)}},simplify=T))


# Loading functions -------------------------------------------------------

source("./real_data/R/run.simulation.R")
source("./real_data/R/Simulation_script.R")

# Setting parameters ------------------------------------------------------

alphas <- c(0.05, 0.1, 0.2, 0.5, 1, 2) 
percentage <- c(0.05,0.1,0.3,0.5,0.7,0.9)
traits <- read.csv("./real_data/data/analysis_data.csv", header = TRUE)
traits[, 2:3] <- log(traits[, 2:3])
tree <- read.newick("./real_data/data/analysis_tree.tree")
dat <- traits[, 2:3]
row.names(dat) <- traits[, 1]
dat <- dat[match(tree$tip.label, row.names(dat)), ]
phylo_sig <- phylolm(Final.Brain.Weight..g. ~ 1, data = dat, phy = tree, model = "OUfixedRoot") 

res <- sim.real(traits, tree = tree, alphas = alphas[1], percentage = percentage[1], method = "dist", correl="rand")

test_that("size of sim.2 results", {
  
  expect_equal(length(res), 2)
  
}
)

test_that("size of sim.2 table and class", {
  
  expect_equal(length(res$statistics), 11)
  expect_equal(class(res$statistics), "numeric")
  expect_equal(dim(res$traits)[2], 4)
  
}
)

