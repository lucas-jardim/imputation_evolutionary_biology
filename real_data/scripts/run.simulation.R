run_simulation <- function(traits, tree, alphas, percentage, output, cl, missing_proc){



if(missing_proc == "rand"){ 
    
clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.real")

dist.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine = "cbind") %:% foreach(1:100,.packages=pacotes,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "dist", correl="rand")

listwise.resul.rand<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "listwise", correl="rand")

mice.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", correl="rand")

mice.phy.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", phylo = T, correl="rand")

PEM.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", correl="rand")

PEM.notrait.resul.rand <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", trait.PEM = F, correl="rand")

obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
stopCluster(clust)

  }
if(missing_proc == "phylo"){
    
clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.real")
dist.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "dist", correl="phylo")

listwise.resul.phylo<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "listwise", correl="phylo")

mice.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", correl="phylo")

mice.phy.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", phylo = T, correl="phylo")

PEM.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", correl="phylo")

PEM.notrait.resul.phylo <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", trait.PEM = F, correl="phylo")

stopCluster(clust)
obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
  }
  
if(missing_proc == "trait"){

clust <- makeCluster(cl)
registerDoParallel(clust)
clusterExport(clust, "sim.real")
dist.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "dist", correl="trait")

listwise.resul.trait<- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "listwise", correl="trait")

mice.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", correl="trait")

mice.phy.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "mice", phylo = T, correl="trait")

PEM.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", correl="trait")

PEM.notrait.resul.trait <- foreach(i=alphas,.packages= pacotes) %:% foreach(j=percentage,.packages=pacotes,.combine="cbind") %:% foreach(1:100,.errorhandling="pass") %dopar% sim.real(traits, tree = tree, alphas = i, percentage = j, method = "PEM", trait.PEM = F, correl="trait")

stopCluster(clust)
obj <- grep("resul", ls(all.names = TRUE), value = TRUE)

save(list = obj, file = output, envir = environment())
  }
  
}
