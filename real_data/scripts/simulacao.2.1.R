
pacotes <- c("mice","geiger","ape","phytools","phylolm","MPSEM","doParallel","letsR","vegan")
invisible(sapply(pacotes,
                 FUN= function(x) {if(!require(x,character.only=T)){install.packages(x,dep=T)
                                                                        require(x)}},simplify=T))

source("./real_data/scripts/run.simulation.R")
source("./real_data/scripts/Simulation_script.R")


alphas <- c(0.05, 0.1, 0.2, 0.5, 1, 2) 
percentage <- c(0.05,0.1,0.3,0.5,0.7,0.9)
cl = 6
traits <- read.csv("./real_data/data/analysis_data.csv", header = TRUE)
traits[, 2:3] <- log(traits[, 2:3])
tree <- read.newick("./real_data/data/analysis_tree.tree")
dat <- traits[, 2:3]
row.names(dat) <- traits[, 1]
dat <- dat[match(tree$tip.label, row.names(dat)), ]
phylo_sig <- phylolm(Final.Brain.Weight..g. ~ 1, data = dat, phy = tree, model = "OUfixedRoot") #mudar para ou estimation no geiger



# run random missing data

dir.create("./real_data/results/0.6/rand", recursive = TRUE)

run_simulation(traits, tree, alphas = alphas, percentage = percentage, output = "./real_data/results/0.6/rand/result.RData", cl = cl, missing_proc = "rand")

# run phylo

dir.create("./real_data/results/0.6/phylo", recursive = TRUE)

run_simulation(traits, tree, alphas = alphas, percentage = percentage, output = "./real_data/results/0.6/phylo/result.RData", cl = cl, missing_proc = "phylo")

# run trait

dir.create("./real_data/results/0.6/trait", recursive = TRUE)

run_simulation(traits, tree, alphas = alphas, percentage = percentage, output = "./real_data/results/0.6/trait/result.RData", cl = cl, missing_proc = "trait")



for(i in c("trait", "phylo", "rand")){
  
  load(paste0("./real_data/results/0.6/", i, "/result.RData"))

}

check.results <- function (resultado) {
                    sapply(1:length(resultado), FUN = 
                          function(z) sapply(1:ncol(resultado[[z]]), FUN = 
                                        function(y) sapply(1:length(resultado[[z]][,y]), FUN = 
                                                     function(x) if(is.null(resultado[[z]][,y][[x]]$statistics)){c(z,y,x)})))
                    
}

#(errors <- do.call("rbind", check.results(dist.resul.phylo)))###encontrar as falhas)


######Construcao da tabela#####################################################################                        
metodos <- c("PEM.notrait","PEM.resul","listwise","dist.resul","mice.phy", "mice.resul")

matriz <- NULL
for(m in 1:length(metodos)){
  metodos.1 <- ls()[grep(metodos[m],ls())]
  cenario <- c("\\.rand","\\.phylo","\\.trait")
  for(k in 1:length(cenario)){
  for (i in 1:length(alphas)){
   for(j in 1:length(percentage)){
    cenario.1 <-  grep(cenario[k],metodos.1)
    matriz1 <- do.call("rbind",lapply(get(metodos.1[cenario.1])[[i]][,j], function(t) return(t$statistics)))
    
    #if(m==3){matriz1 <- cbind(matriz1[,1:2],rep(NA,nrow(matriz1)),matriz1[,3:5])}
    matriz1 <- cbind(matriz1,rep(metodos[m],nrow(matriz1)),rep(cenario[k],nrow(matriz1)),rep(alphas[i],nrow(matriz1)),rep(percentage[j],nrow(matriz1)))
    print(c(dim(matriz), dim(matriz1)))
    matriz <- rbind(matriz,matriz1)
    cat(c(m,k,i,j),"\n")
    }
}
}
}

matriz.1 <- as.data.frame(matriz,stringsAsFactors=F)
matriz.1[,1:11] <- apply(matriz[,-(12:15)],2,FUN=as.numeric)
colnames(matriz.1)[12:ncol(matriz.1)] <- c("methods", "mechanism", "alpha", "percentage")
colnames(matriz.1)[11] <- "NRMSE"
matriz.1[, "methods"] <- gsub("\\.resul", "", matriz.1[, "methods"])
matriz.1[, "methods"] <- gsub("dist", "hot-deck", matriz.1[, "methods"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.rand", "MCAR", matriz.1[, "mechanism"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.phylo", "MAR-PHYLO", matriz.1[, "mechanism"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.trait", "MAR-TRAIT", matriz.1[, "mechanism"])
matriz.1[, "alpha"] <- rep(phylo_sig$optpar, length(matriz.1$alpha))
metodos <- c("PEM.notrait","PEM","listwise","hot-deck","mice.phy", "mice")
cenario <- c("MCAR", "MAR-PHYLO", "MAR-TRAIT")
