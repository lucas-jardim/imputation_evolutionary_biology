library(RColorBrewer)

grafico <- function(alphas,cenario,ylimit,ylimiti,matriz,coluna,cores,space, yleg){
  
  ###########################################################################
  # alphas = A numeric vector with the simulated values of Ornstein-Uhlenbeck's alphas.
  # cenario = A character vector with simulated missing data mechanisms.
  # ylimit = Plot's upper limits of y axis.
  # ylimiti = Plot's lower limits of y axis.
  # matriz = Results matrix.
  # coluna = Matrix column (matriz) to be plotted.
  # cores = A vector of color.
  # space = Steps in axis  values.
  # yleg = A character vector with simulated missing data mechanisms.
  ###########################################################################
  
  
    for(j in 1:length(cenario)){
      plot.new()
      plot.window(xlim=c(0.05,1),ylim=c(ylimiti,ylimit))
      axis(1,at=percentage,cex=1.5)   
      axis(2,pos=c(0.05,0),seq(ylimiti,ylimit,space),cex=0.1)
      text(substitute(paste(alpha,"=",valor,sep=" "),list(valor=round(as.numeric(alphas), 2))),x=0.15,y=ylimit-0.1*ylimit, cex = 1.5)
      
      for(k in 1:length(metodos)){    
        linhas <- which(matriz[, "methods"]== metodos[k]&matriz[,"mechanism"]==cenario[j]&matriz[,"alpha"]==as.character(alphas))
        percentage_groups <- split(matriz[linhas, c(coluna - 1, coluna)], as.factor(matriz[linhas,"percentage"]))
        statis_diff <- lapply(1:length(percentage_groups), function(g) (percentage_groups[[g]][, 2] - percentage_groups[[g]][, 1])/percentage_groups[[g]][, 1])
        names(statis_diff) <- names(percentage_groups)
        quantile_stats <- sapply(statis_diff, quantile, 0.5)
        lines(as.numeric(names(statis_diff)), quantile_stats, col = cores[k])
      }
      lines(seq(0,0.9,0.1),rep(0,length(seq(0,0.9,0.1))),cex=2,col="black",lty=2) 
      if(j==1){legend(0.75,ylimit-0.1*ylimit,metodos,fill=cores,cex=1,bty="n")}
      mtext(yleg, side = 2, outer = TRUE, line = 1)
      mtext("Percentage of missing data", side = 1, outer = TRUE, line = 2)
    }
  
  mtext("MCAR",outer=T,at=c(0.17,4))
  mtext("MAR-PHYLO",outer=T,at=c(0.5,4))
  mtext("MAR-TRAIT",outer=T,at=c(0.85,4))
  par(mfrow=c(1,1))
}
cores <- RColorBrewer::brewer.pal(6, "Dark2")
yleg <- c("Blomberg's K error", "Moran's I error", "Mean error", "Variance error", "Regression Coefficient error")
dir.create("./real_data/graphics")
count <- 1
for(i in seq(2, 10, 2)){
tiff(filename=paste0("./real_data/graphics/",colnames(matriz.1)[i],".tif"),res=300,width=30,height=30,units="cm")
par(mfrow=c(1,3),mar=c(2,2,1,0),oma=c(3,3,2,0),bty="l")
grafico(as.character(phylo_sig$optpar),cenario,ylimit=2, ylimiti=-1,matriz.1,coluna=i,cores,space=0.25, yleg = yleg[count])
count <- count + 1
dev.off()
}


# Regression tree -----------------------------------------------------------

for(i in c("moran", "beta", "mean", "var", "k")){
  tiff(filename=paste0("./real_data/graphics/", i, "_RT.tif"),res=300,unit="cm",width=30,height=30)
  rpart.plot(get(i), type = 5, cex = 0.9, shadow.col = "grey", branch.lwd = 2)
  dev.off()
}



# Phylogenetic signal vs imputation errors --------------------------------

tiff(filename=paste0("./real_data/graphics/imputation_error.tif"),res=300,unit="cm",width=30,height=30)
par(mfrow = c(1, 2))

plot(x = matriz.1[, "NRMSE"], y = abs((matriz.1[, 2] - matriz.1[, 1])/matriz.1[, 1]), xlab = "Imputation error", ylab= "Blomberg's K error" ,pch = 19) 
text(1.3, 8, "A")
plot(x = matriz.1[, "NRMSE"], y = abs((matriz.1[, 4] - matriz.1[, 3])/matriz.1[, 3]), xlab = "Imputation error", ylab= "Moran's I error", pch = 19) 
text(1.3, 1.78, "B")
dev.off()


