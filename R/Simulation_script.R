sim.2 <- function(sp,alphas,percentage,method,r,m=10,phylo=F,trait.PEM=T,correl, sim_phylo = pbtree(n = sp)){

###########################################################################  
# sp = number of species
# alphas = Ornstein-Uhlenbeck model's parameter.
# percentage = percentabe of missing data.  
# method = method to dealing missing data. 
# r = Correlation between traits.
# m = Amount of multiple imputation datasets.
# phylo = Boolean argument. If FALSE phylogenetic eigenvectors are not included in mice.
# trait.PEM = Boolean argument. If TRUE an explanatory variable is included in PEM.
# correl = Mechanism of missing data. The options are "rand", "phylo", "trait".
# sim_phylo = Simulated phylogeny. The default is phytools::pbtree(n = sp). 
###########################################################################  
  
# Simulating phylogeny and traits -----------------------------------------

  tree <- sim_phylo
  tree$node.label <- paste("N", 1:tree$Nnode, sep = "")
  trait <- rTraitCont(phy=tree,model="OU",sigma = 1, alpha =alphas, theta = 0) # criar traits nas filogenias
  trait2 <-if(correl!="trait"){ r*trait+sqrt(1-r^2)*fastBM(tree,sig2=mean(pic(trait,tree)^2))}
  
# Generating missing values -----------------------------------------------

  if(correl=="rand"){
    # excluding values by chance
    
    trait1 <- trait
    trait1[sample(1:length(trait1),percentage*length(trait1))] <- NA 
    na.sp <- is.na(trait1)
  }
  if(correl=="phylo"){
    # Phylogenetically correlated missing values
    
    distancia.phy <- cophenetic(tree)
    distancia.phy <- distancia.phy[sample(1:nrow(distancia.phy),1),]
    trait1<- trait
    trait1[order(distancia.phy)[1:(percentage*length(trait1))]] <- NA    
    na.sp <- is.na(trait1)   
  }
  
  if(correl=="trait"){
    # Missing data correlated to a trait
    
    COV <- matrix(0,length(tree$tip.label),length(tree$tip.label))
    diag(COV) <- diag(vcv(tree))
    trait2 <-  r*trait+sqrt(1-r^2)*mvrnorm(1,mu=rep(0,length(tree$tip.label)),Sigma=COV)
    trait1 <- trait    
    trait1[order(trait2)[1:(percentage*length(trait1))]] <- NA    
    na.sp <- is.na(trait1)
    
  }

# Imputation --------------------------------------------------------------

  if(method=="dist"){
    # Imputing by phylogenetic distance
    
    distance <- cophenetic(tree)
    diag(distance) <- NA
    missing_dist <- distance[names(trait1[na.sp]), !na.sp]
    closest_sp <- apply(missing_dist, 1, function(i) names(which(i == min(i, na.rm = TRUE))))
    source_sp <- sapply(closest_sp, function(i) sample(i, 1))
    trait1[na.sp] <- trait1[source_sp] # imputed data 
    ktrait <- phylosig(tree=tree,x=trait) # phylogenetic signal of complete data 
    ktrait1 <- phylosig(tree=tree,x=trait1) # phylogenetic signal of missing data
    distancia.phy <- cophenetic(tree)
    sturges <- 1+3.3*log((length(trait)*(length(trait)-1))/2)
    moran <- lets.correl(trait,distancia.phy,sturges,equidistant=T,plot=F)
    moran <- moran[1, 1]
    moran1 <- lets.correl(trait1,distancia.phy,sturges,equidistant=T,plot=F)
    moran1 <- moran1[1, 1]
    sd.mean <- sqrt(mean((trait[na.sp]-trait1[na.sp])^2)/(max(trait)-min(trait)))
    mean.trait<- mean(trait)
    mean.trait1 <- mean(trait1)
    var.trait <- var(trait)
    var.trait1 <- var(trait1)
    reg <- lm(trait ~ trait2)
    beta.trait <- coef(reg)[2]
    reg1 <- lm(trait1 ~ trait2)
    beta.trait1 <- coef(reg1)[2]    
    statistics <- c(ktrait, ktrait1, moran, moran1, mean.trait, mean.trait1, 
                    var.trait, var.trait1, beta.trait, beta.trait1,sd.mean)
    names(statistics) <- c("ktrait", "ktrait1", "moran", "moran1", 
                           "mean.trait", "mean.trait1", "var.trait", "var.trait1", 
                           "beta.trait", "beta.trait1", "sd.mean")
    
  }
  
  if(method =="mice"){
    # Multiple imputation
    
    matriz <- cbind(trait,trait1)# matrix with complete and excluded data
    matriz.mice <- cbind(trait1,trait2)
    if(phylo==F){mice1 <-mice(matriz.mice,m)}else{
      pruned <- drop.tip(tree, names(trait1[na.sp]))
      graphs <-  Phylo2DirectedGraph(pruned)
      fit <- PEM.fitSimple(y = trait1[na.sp==F], x=trait2[na.sp==F], w = graphs,lower=0,upper=1)
      PEM.t <- PEM.build(graphs, d="distance",sp="species",a=unlist(fit$a), psi=unlist(fit$psi))                                                            
      eigens.t <- lmforwardsequentialAICc(y = trait1[na.sp==F],x=trait2[na.sp==F],object=PEM.t)
      loc <- getGraphLocations(tree,names(trait1[na.sp]))
      scores <- Locations2PEMscores(fit, loc)
      PEM.t <- rbind(PEM.t$u,scores$scores)
      PEM.t <- PEM.t[row.names(matriz.mice),]
      matriz.mice <- cbind(matriz.mice,PEM.t[,attr(eigens.t$terms,"term.labels")[-1]])
      mice1 <- mice(matriz.mice,m) 
    }
    mice1.complete <- sapply(1:mice1$m, FUN = function(y) complete(mice1,y)[,1])
    ktrait <- phylosig(tree=tree,x=trait)
    ktrait1 <- mean(apply(mice1.complete,2,FUN = phylosig,tree=tree))
    distancia.phy <- cophenetic(tree)
    sturges <- 1+3.3*log((length(trait)*(length(trait)-1))/2)
    moran <- lets.correl(trait,distancia.phy,sturges,equidistant=T,plot=F)
    moran <- moran[1, 1]
    moran1 <- sapply(1:ncol(mice1.complete),FUN= function(x)lets.correl(mice1.complete[,x],y=distancia.phy,z=sturges,equidistant=T,plot=F),simplify=F)
    moran1 <- mean(sapply(moran1, function(x) x[1,1]))
    sd.mean <- sqrt(mean((trait[na.sp]-rowMeans(mice1.complete[na.sp,]))^2)/(max(trait)-min(trait)))###NRMSE ,permite calcular a eficiencia da imputacao como Penone
    mean.trait <- mean(trait)
    mean.trait1 <- mean(apply(mice1.complete, 2, mean))
    var.trait <- var(trait)
    var.trait1 <- mean(apply(mice1.complete, 2, var))
    reg <- lm(trait ~ trait2)
    beta.trait <- coef(reg)[2]
    reg1 <- apply(mice1.complete, 2, FUN = function(x) lm(x ~ trait2))
    beta.trait1 <- mean(sapply(reg1, function(x) coef(x)[2]))
    statistics <- c(ktrait, ktrait1, moran, moran1,mean.trait, mean.trait1, 
                    var.trait, var.trait1, beta.trait, beta.trait1, sd.mean)
    names(statistics) <- c("ktrait", "ktrait1", "moran", "moran1",  
                           "mean.trait", "mean.trait1", "var.trait", "var.trait1", 
                           "beta.trait", "beta.trait1", "sd.mean")
    
  }
  
  if(method == "listwise"){
    # Ignoring missing data
    
    tree1 <- drop.tip(tree,names(trait1)[na.sp])
    trait1.s <- trait1[!na.sp]
    ktrait <- phylosig(tree=tree,x=trait)
    ktrait1 <- phylosig(tree=tree,x=trait1.s)
    distancia.phy <- cophenetic(tree)
    distancia.phy1 <- cophenetic(tree1)
    sturges <- 1+3.3*log((length(trait)*(length(trait)-1))/2)
    moran <- lets.correl(trait,distancia.phy,sturges,equidistant=T,plot=F)
    moran <- moran[1,1]
    moran1 <- lets.correl(trait1.s,distancia.phy1,sturges,equidistant=T,plot=F)
    moran1 <- moran1[1,1]
    sd.mean <- NA
    mean.trait <- mean(trait)
    mean.trait1 <- mean(trait1.s)
    var.trait <- var(trait)
    var.trait1 <- var(trait1.s)
    reg <- lm(trait ~ trait2)
    beta.trait <- coef(reg)[2]
    trait2.1 <- trait2[!na.sp]
    reg1 <- lm(trait1.s ~ trait2.1)
    beta.trait1 <- coef(reg1)[2]    
    statistics <- c(ktrait, ktrait1, moran, moran1, mean.trait, mean.trait1, 
                    var.trait, var.trait1, beta.trait, beta.trait1,sd.mean)
    names(statistics) <- c("ktrait", "ktrait1", "moran", "moran1",  
                           "mean.trait", "mean.trait1", "var.trait", "var.trait1", 
                           "beta.trait", "beta.trait1", "sd.mean")
  }
  
   if(method == "PEM"){
    # Phylogenetic eigenvector mapping
     
    matriz <- cbind(trait,trait1)
    pruned <- drop.tip(tree, names(trait1[na.sp]))
    graphs <-  Phylo2DirectedGraph(pruned)
    if(trait.PEM==F){
      fit <- PEM.fitSimple(y = trait1[na.sp==F],x=NULL, w = graphs,lower=0,upper=1)
      PEM.t <- PEM.build(graphs, d="distance",sp="species",a=unlist(fit$a), psi=unlist(fit$psi))
      eigens.t <- lmforwardsequentialAICc(y = trait1[na.sp==F],object=PEM.t)
      loc <- getGraphLocations(tree,names(trait1[na.sp])) 
      trait1[na.sp] <- predict(fit,loc,eigens.t)
      
    }else{
      fit <- PEM.fitSimple(y = trait1[na.sp==F], x=trait2[na.sp==F], w = graphs,lower=0,upper=1)
      PEM.t <- PEM.build(graphs, d="distance",sp="species",a=unlist(fit$a), psi=unlist(fit$psi))
      eigens.t <- lmforwardsequentialAICc(y = trait1[na.sp==F], x = trait2[na.sp==F], object=PEM.t)
      loc <- getGraphLocations(tree,names(trait1[na.sp]))
      newdata <- as.matrix(trait2[na.sp])
      colnames(newdata) <- "x"
      trait1[na.sp] <- predict(object = fit, targets = loc, lmobject = eigens.t, newdata= newdata)
    } 
    
    ktrait <- phylosig(tree=tree,x=trait)
    ktrait1 <- phylosig(tree=tree,x=trait1)
    distancia.phy <- cophenetic(tree)
    sturges <- 1+3.3*log((length(trait)*(length(trait)-1))/2)
    moran <- lets.correl(trait,distancia.phy,sturges,equidistant=T,plot=F)
    moran <- moran[1, 1]
    moran1 <- lets.correl(trait1,distancia.phy,sturges,equidistant=T,plot=F)
    moran1 <- moran1[1, 1]
    sd.mean <- sqrt(mean((trait[na.sp]-trait1[na.sp])^2)/(max(trait)-min(trait)))
    mean.trait<- mean(trait)
    mean.trait1 <- mean(trait1)
    var.trait <- var(trait)
    var.trait1 <- var(trait1)
    reg <- lm(trait ~ trait2)
    beta.trait <- coef(reg)[2]
    reg1 <- lm(trait1 ~ trait2)
    beta.trait1 <- coef(reg1)[2]    
    statistics <- c(ktrait, ktrait1, moran, moran1, mean.trait, mean.trait1, 
                    var.trait, var.trait1, beta.trait, beta.trait1, sd.mean)
    names(statistics) <- c("ktrait", "ktrait1", "moran", "moran1",  
                           "mean.trait", "mean.trait1", "var.trait", "var.trait1", 
                           "beta.trait", "beta.trait1", "sd.mean")
    
  }   
  
  return(list( statistics = statistics, traits = cbind(trait, trait1, trait2, na.sp)))
} 