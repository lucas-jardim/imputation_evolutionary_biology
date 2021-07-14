library(phytools)
library(letsR)
library(testthat)

source("./R/Simulation_script.R")

alphas= 1.5
percentage=0.5
r=0.6
m=10
correl="rand"
sp=200
phylo=T

res <- sim.2(sp = sp, 
          alphas = alphas,
          percentage = percentage, 
          method = "dist",
          r = r,
          m = m,
          phylo = FALSE,
          trait.PEM = TRUE,
          correl="rand"
          )


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




