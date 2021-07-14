run_simulation <- function(alphas, percentage, r, sp, output, cl, missing_proc){

###########################################################################
# alphas = Ornstein-Uhlenbeck's alpha.
# percentage = Percentage of missing data.
# r = Trait's correlation.
# sp = Number of species.
# output = Path to save output.
# cl = Number of cluster in parallel computing.
# missing_proc = Mechanism of missing data. The options are "rand", "phylo" 
#                and " trait".
###########################################################################

if(missing_proc == "rand"){ 
    

# Setting clusters --------------------------------------------------------

clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.2")

# Running simulations -----------------------------------------------------

dist.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine = "cbind") %:% foreach(1:100,.packages=pacotes,.errorhandling="pass") %dopar% sim.2(sp,i,j,"dist",correl="rand",r)

listwise.resul.rand<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp , i,j,"listwise",correl="rand",r)

mice.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",correl="rand",r)

mice.phy.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",phylo=T,correl="rand",r)

PEM.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"PEM",correl="rand",r)

PEM.notrait.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar%sim.2(sp,i,j,"PEM",trait.PEM=F,correl="rand",r)

obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
stopCluster(clust)

  }
if(missing_proc == "phylo"){
    
  # Setting clusters --------------------------------------------------------
clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.2")
dist.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"dist",correl="phylo",r)

# Running simulations -----------------------------------------------------

listwise.resul.phylo<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"listwise",correl="phylo",r)

mice.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",correl="phylo",r)

mice.phy.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",phylo=T,correl="phylo",r)

PEM.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"PEM",correl="phylo",r)

PEM.notrait.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar%sim.2(sp,i,j,"PEM",trait.PEM=F,correl="phylo",r)

stopCluster(clust)
obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
  }
  
if(missing_proc == "trait"){

# Setting clusters --------------------------------------------------------

clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.2")


# Running simulations -----------------------------------------------------

dist.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"dist",correl="trait",r)

listwise.resul.trait<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"listwise",correl="trait",r)

mice.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",correl="trait",r)

mice.phy.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"mice",phylo=T,correl="trait",r)

PEM.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.2(sp,i,j,"PEM",correl="trait",r)

PEM.notrait.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar%sim.2(sp,i,j,"PEM",trait.PEM=F,correl="trait",r)

stopCluster(clust)

obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
  }
  
}
