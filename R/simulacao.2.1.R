pacotes <- c("mice","geiger","ape","phytools","MPSEM","doParallel","letsR","vegan")


# Importing packages ------------------------------------------------------

invisible(sapply(pacotes,
                 FUN= function(x) {if(!require(x,character.only=T)){install.packages(x,dep=T)
                                                                        require(x)}},simplify=T))

# Importing functions -----------------------------------------------------

source("./R/run.simulation.R")
source("./R/Simulation_script.R")

#alphas= 1.5
#percentage=0.5
#r=0.6
#m=10
#correl="rand"
#sp=200
#phylo=T
#phylo = pbtree( n = sp)


# Setting parameters ------------------------------------------------------

alphas <- c(0.05, 0.1, 0.2, 0.5, 1, 2)
percentage <- c(0.05,0.1,0.3,0.5,0.7,0.9)
r <- 0.6 
sp <- 200

# Run random missing data -------------------------------------------------

dir.create("./Results/0.6/rand", recursive = TRUE)

run_simulation(alphas = alphas, percentage = percentage, r= r, sp = sp, output = "./Results/0.6/rand/result.RData", cl = 10, missing_proc = "rand")

# run phylo

# Run phylogenetically correlated missing data ----------------------------

dir.create("./Results/0.6/phylo", recursive = TRUE)

run_simulation(alphas = alphas, percentage = percentage, r= r, sp = sp, output = "./Results/0.6/phylo/result.RData", cl = 10, missing_proc = "phylo")


# Run missing data correlated to a trait ----------------------------------

dir.create("./Results/0.6/trait", recursive = TRUE)

run_simulation(alphas = alphas, percentage = percentage, r= r, sp = sp, output = "./Results/0.6/trait/result.RData", cl = 10, missing_proc = "trait")


# Importing simulations results -------------------------------------------

for(i in c("trait", "phylo", "rand")){
  
  load(paste0("./Results/0.6/", i, "/result.RData"))

}


# Checking simulation errors ----------------------------------------------


check.results <- function (resultado) {
                    sapply(1:length(resultado), FUN = 
                          function(z) sapply(1:ncol(resultado[[z]]), FUN = 
                                        function(y) sapply(1:length(resultado[[z]][,y]), FUN = 
                                                     function(x) if(is.null(resultado[[z]][,y][[x]]$statistics)){c(z,y,x)})))
                    
}

#(errors <- do.call("rbind", check.results(mice.resul.trait)))###encontrar as falhas)
#lapply(1:nrow(errors), function(x) mean.resul.trait[[errors[x,1]]][,errors[x,2]][[errors[x,3]]])

#PEM.resul.trait[[errors[,1]]][,errors[,2]][[errors[,3]]] <- sim.2(sp,alphas[2],percentage[1],"PEM",correl="trait",r)


# Creating result's table -------------------------------------------------

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


# Cleaning result's table -------------------------------------------------

matriz.1 <- as.data.frame(matriz,stringsAsFactors=F)
matriz.1[,1:11] <- apply(matriz[,-(12:15)],2,FUN=as.numeric)
colnames(matriz.1)[12:ncol(matriz.1)] <- c("methods", "mechanism", "alpha", "percentage")
colnames(matriz.1)[11] <- "NRMSE"
matriz.1[, "methods"] <- gsub("\\.resul", "", matriz.1[, "methods"])
matriz.1[, "methods"] <- gsub("dist", "hot-deck", matriz.1[, "methods"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.rand", "MCAR", matriz.1[, "mechanism"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.phylo", "MAR-PHYLO", matriz.1[, "mechanism"])
matriz.1[, "mechanism"] <- gsub("\\\\\\.trait", "MAR-TRAIT", matriz.1[, "mechanism"])

metodos <- c("PEM.notrait","PEM","listwise","hot-deck","mice.phy", "mice")
cenario <- c("MCAR", "MAR-PHYLO", "MAR-TRAIT")
