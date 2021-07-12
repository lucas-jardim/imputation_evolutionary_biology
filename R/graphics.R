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
  
  for(i in c(1, 4, 6)){
    for(j in 1:length(cenario)){
      plot.new()
      plot.window(xlim=c(0.05,1),ylim=c(ylimiti,ylimit))
      axis(1,at=percentage,cex=1.5)   
      axis(2,pos=c(0.05,0),seq(ylimiti,ylimit,space),cex=0.1)
      text(substitute(paste(alpha,"=",valor,sep=" "),list(valor=alphas[i])),x=0.15,y=ylimit-0.1*ylimit)
      
      for(k in 1:length(metodos)){    
        linhas <- which(matriz[, "methods"]== metodos[k]&matriz[,"mechanism"]==cenario[j]&matriz[,"alpha"]==as.character(alphas[i]))
        percentage_groups <- split(matriz[linhas, c(coluna - 1, coluna)], as.factor(matriz[linhas,"percentage"]))
        statis_diff <- lapply(1:length(percentage_groups), function(g) (percentage_groups[[g]][, 2] - percentage_groups[[g]][, 1])/percentage_groups[[g]][, 1])
        names(statis_diff) <- names(percentage_groups)
        quantile_stats <- sapply(statis_diff, quantile, 0.5)
        lines(as.numeric(names(statis_diff)), quantile_stats, col = cores[k])
      }
      lines(seq(0,0.9,0.1),rep(0,length(seq(0,0.9,0.1))),cex=2,col="black",lty=2) 
      if(i==1&j==1){legend(0.5,ylimit-0.1*ylimit,metodos,fill=cores,cex=0.7,bty="n")}
      mtext(yleg, side = 2, outer = TRUE, line = 1)
      mtext("Percentage of missing data", side = 1, outer = TRUE, line = 2)
    }
  }
  mtext("MCAR",outer=T,at=c(0.17,4))
  mtext("MAR-PHYLO",outer=T,at=c(0.5,4))
  mtext("MAR-TRAIT",outer=T,at=c(0.85,4))
  par(mfrow=c(1,1))
}

# Plotting simulation results ---------------------------------------------

cores <- RColorBrewer::brewer.pal(6, "Dark2")
yleg <- c("Blomberg's K error", "Moran's I error", "Mean error", "Variance error", "Regression Coefficient error")
dir.create("./graphics")
count <- 1

for(i in seq(2, 10, 2)){
    tiff(filename=paste0("./graphics/",colnames(matriz.1)[i],".tif"),res=300,width=30,height=30,units="cm")
    par(mfrow=c(3,3),mar=c(2,2,1,0),oma=c(3,3,2,0),bty="l")
    grafico(alphas,cenario,ylimit=5, ylimiti=-1,matriz.1,coluna=i,cores,space=0.25, yleg = yleg[count])
    count <- count + 1
    dev.off()
}

# Plotting regression tree ------------------------------------------------

for(i in c("moran", "beta", "mean", "var")){
   tiff(filename=paste0("./graphics/", i, "_RT.tif"),res=300,unit="cm",width=30,height=30)
   rpart.plot(get(i), type = 5, cex = 0.9, shadow.col = "grey", branch.lwd = 2)
   dev.off()
}

