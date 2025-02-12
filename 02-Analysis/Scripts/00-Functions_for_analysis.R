# beneat marine 
# 14/10/24 
# functions used by the annex scripts

library(knitr)
library(kableExtra)
library(archetypes)
library(ggtern)
library(tidyverse)
library(archetypes)
library(stringr)
library(dplyr)
library(vegan) #decostand
library(boot)
library(viridis)
library(ggtern)
library(Ckmeans.1d.dp)
library(FactoMineR)
library(ggrepel)
library(ggfortify)
library(RColorBrewer)
library(factoextra)
library(ape)
library(phytools)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rphylopic)
library(ggforce)
library(ggtree)
library(phylopath)
library(scales)
library(plotROC)
library(pROC)
library(rstatix)
library(ggpubr)

suppressWarnings(RNGversion("3.5.0"))
set.seed(1986)
mycol = function(x,alpha){
  y = col2rgb(x)
  rgb( y[1,1]/255, y[2,1]/255, y[3,1]/255, alpha=alpha)
}


addhabitatcolumn <- function(dataphylo){
  hhabtot <- apply(dataphylo[, c("habitatpelagic", "habitatdemersal", "habitatbenthopelagic")], 1, which.max)
  dataphylo$hhabtot <- hhabtot
  return(dataphylo)
}


ratioTOtrait <- function(ratio, operator, IS_OPERATOR_Loo, IS_LOG_M){
  
  if (IS_LOG_M){
    if (IS_OPERATOR_Loo)  {
      trait <- ratio + operator
    }
    else {
      operator[which(operator<0)] <-0.01
      trait <- ratio - operator
    }
  }
  else {
    if (IS_OPERATOR_Loo)  {
      trait <- ratio + operator
    }
    else {
      operator[which(operator<0)] <-0.01
      trait <- log(exp(ratio) / operator)
    }
  }
  return(trait)
}



runAA <- function(dataplot, traits, kmax){
  beta_iv1 = dataplot[, which(colnames(dataplot) %in%  traits)]
  repAA<-stepArchetypes(data=beta_iv1, k=1:kmax, verbose=TRUE, nrep=30)
  AA<-bestModel(repAA[[min(kmax, 3)]])
  return(AA)
}



plotggtern <- function(dataphylo, AA, vecAA_eq, bins){
  cmstd <- dataphylo$c_m
  df <- cbind(AA$alphas, cmstd)
  colnames(df) <- c(vecAA_eq, "d")
  dr <- df
  df <- as.data.frame(df)
  plot<- ggtern::ggtern(df, aes(x=a, y=b, z=c))+ 
    ggtern::geom_tri_tern(bins = bins, fun = median, aes(value = d, fill = after_stat(stat))) +  
    stat_tri_tern(bins = bins, fun = median,
                  geom = 'text',
                  aes(value = d, label = sprintf("%.2f", after_stat(stat))),
                  size = 4, color = 'black', centroid = TRUE, fontface = "bold") +
    
    # Use scale_fill_gradientn() for custom color gradient
    scale_fill_gradientn(colors = c("darkorchid4", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027"), 
                         name = "Standardized \nmass-specific \nmaintenance")
  return(plot)
}




dataacp_add_colorvector <- function(dataphylo, kclusters, dataacp){
  
  # clusters with k-means
  x <- dataphylo$c_m
  result <- Ckmeans.1d.dp(x, k=1:kclusters)
  k <- max(result$cluster)
  centers <- result$centers
  cluster_centroid <- centers[result$cluster]
  
  # extract clusters to use as color vector
  colACP <- cluster_centroid
  dataacp <- cbind(dataacp, colACP)
  quantiles = as.vector(quantile(dataphylo$c_m, seq(0,1,0.1), na.rm=T)) 
  dataacp$colquantiles <- cut(dataphylo$c_m, breaks = unique(quantiles), include.lowest = TRUE)
  
  # color or number per Family or class
  if (dim(table(dataphylo$Class))>2){
    genecolour <- rep('NA', nrow(df))
    genecolour[dataphylo$Family == 'Gobiidae'] <- 'brown3'
    genecolour[dataphylo$Family == 'Clupeidae'] <- 'forestgreen'
    genecolour[dataphylo$Class == 'Elasmobranchii'] <- 'royalblue'
    dataacp$genecolour <- genecolour
    
    dataacp$colClass <- as.numeric(as.factor(dataacp$Class))
    dataacp$colClass[which(dataacp$colClass==1)] <- c(8)
    dataacp$colClass[which(dataacp$colClass==2)] <- c(25)
    dataacp$colClass[which(dataacp$colClass==3)] <- c(19)
  } 
  
  if (dim(table(dataphylo$Class))>2){
    #color examples 
    famOPP <- c("Gobiidae", "Engraulidae", "Atherinopsidae", "Fundulidae", "Gambusia")
    famELA <- c("Ariidae")
    famPER <- c("Gadidae", "Moridae", "Melanonidae", "Euclichthyidae", "Rachycentridae", "Molidae", "")
    
    dataacp$colFAM <- ifelse(dataphylo$Family==famELA, c("equi"), "else")
    dataacp$colFAM[dataphylo$Family==famOPP] <- "opp"
    dataacp$colFAM[dataphylo$Family==famPER] <- "per"
  } 
  return(list(dataacp, result$centers))
}



runPCA <- function(dataplot, traits){
  beta_iv1 = dataplot[, which(colnames(dataplot) %in%  traits)]
  PCA <- prcomp(beta_iv1) 
  return(PCA)
}


runpPCA <- function(dataplot, traits, tree){
  beta_iv1 = dataplot[, which(colnames(dataplot) %in%  traits)]
  PCA <- phyl.pca(tree, beta_iv1, method="BM", mode="cov")
  return(PCA)
}


preparedataforplot <- function(numbPCA1, numbPCA2, dataacp, AA, PCA){
  dataacp <- dataacp[[1]]
  dataacp$PC2 <- PCA$x[, numbPCA2] # indexing the first column
  dataacp$PC1 <- PCA$x[, numbPCA1]  # indexing the second column
  rotation_data <- PCA$rotation[, c(numbPCA1, numbPCA2)]
  colnames(rotation_data) <- c(paste0("PC", AXESTOREPRESENT[1]), paste0("PC", AXESTOREPRESENT[2]))
  archetypepointsPCA <- scale(AA$archetypes, PCA$center, PCA$scale) %*% PCA$rotation
  eigenval <- get_eigenvalue(PCA)[c(numbPCA1, numbPCA2),"variance.percent"]
  return(list(rotate=rotation_data, AArotate=archetypepointsPCA, dataacp, eigenval))
}

preparedataforplotp <- function(numbPCA1, numbPCA2, dataacp, AA, PCA){
  dataacp <- as.data.frame(dataacp) #[[1]]
  dataacp$PC2 <- PCA$S[, numbPCA2] # indexing the first column
  dataacp$PC1 <- PCA$S[, numbPCA1]  # indexing the second column
  rotation_data <- PCA$L[, c(numbPCA1, numbPCA2)]
  colnames(rotation_data) <- c("PC1", "PC2")
  archetypepointsPCA <- scale(AA$archetypes, PCA$a, FALSE) %*% PCA$L
  eigen <- colSums(PCA$Evec)
  eigenval <- c(abs(colSums(PCA$Evec))/sum(abs(eigen)))[c(numbPCA1, numbPCA2)]
  return(list(rotate=rotation_data, AArotate=archetypepointsPCA, dataacp, eigenval))
}


plotPCAclean <- function(veccol, vecAA, PCA){
  
  palette = "viridis"
  
  plotPCA <- ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_point(data = dataacpPLOT, aes(x = PC1, y = PC2, color = as.factor(colACP)), size=1.5)+ 
    geom_segment(data = rotation,
                 aes(x = 0, y = 0, xend = PC1*pozarrow, yend = PC2*pozarrow),
                 arrow = arrow(length = unit(0.25, "mm")),
                 color = "white",
                 linewidth=1)+
    geom_text_repel(data = rotation,
                    aes(x = PC1*pozarrow+0.22*PC1, y = PC2*pozarrow, label = rownames(rotation)),
                    colour = "black", fontface = c("bold"),hjust=1, size=4.8)+
    geom_text_repel(data=rotation,
              aes(x = PC1*pozarrow+0.22*PC1, y = PC2*pozarrow, label = rownames(rotation)),
              colour = "white", fontface = c("bold"),hjust=1, size=4.5)+
    geom_text_repel() +
    scale_colour_viridis_d(name="c_m clusters", option=palette, direction = -1)+
    geom_point(size=6.2, aes(matrixAAinPCA[,numbPCA1]+0.02, matrixAAinPCA[,numbPCA2]), 
               colour="white")+
    geom_point(size=6, aes(matrixAAinPCA[,numbPCA1], matrixAAinPCA[,numbPCA2]), 
               colour=veccol)+
    geom_text(aes(matrixAAinPCA[,numbPCA1]-0.02-1, c(matrixAAinPCA[,numbPCA2]-0.2-0.01)), 
              colour=c("black"),
              label = vecAA, fontface = c("bold"), hjust=-0.25, size=5.2)+
    geom_text(aes(matrixAAinPCA[,numbPCA1]-1, c(matrixAAinPCA[,numbPCA2]-0.2-0.01)), 
              colour=veccol,
              label = vecAA, fontface = c("bold"), hjust=-0.25, size=5)+
    labs(color="c_m quantiles clustering")+
    xlab(label = paste0("PC", numbPCA1, " ", round(eigenval[1],digits = 2) , "% var"))+
    ylab(label = paste0("PC", numbPCA2, " ", round(eigenval[2],digits = 2) , "% var"))+ 
    theme(
      panel.background = element_rect(fill = "azure3",
                                      colour = "azure3",
                                      linewidth = 0.5, linetype = "solid"),
       panel.grid.minor = element_blank())
  return(list(plotPCA))
}



plotPCA_knownc_m <- function(veccol, vecAA, PCA){
  
  palette = "viridis"
  
  plotPCA <- ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_point(data = dataacpPLOT, aes(x = PC1, y = PC2, color = ifelse(KNOWN, as.factor(colACP), "grey")), 
               size=ifelse(dataacpPLOT$KNOWN, 2, 1.5), 
               alpha = ifelse(dataacpPLOT$KNOWN, 1, 0.05))+ 
    geom_segment(data = rotation,
                 aes(x = 0, y = 0, xend = PC1*pozarrow, yend = PC2*pozarrow),
                 arrow = arrow(length = unit(0.25, "mm")),
                 color = "white",
                 linewidth=1)+
    geom_text(data=rotation,
              aes(x = PC1*pozarrow+0.22, y = PC2*pozarrow, label = rownames(rotation)),
              colour = "white", fontface = c("bold"),hjust=1, size=4)+
    geom_text(data = rotation,
              aes(x = PC1*pozarrow+0.25, y = PC2*pozarrow, label = rownames(rotation)),
              colour = "black", fontface = c("bold"),hjust=1, size=4)+
    geom_text_repel() +
    scale_colour_viridis_d(name="c_m clusters", option=palette, direction = -1)+
    geom_point(size=6, aes(matrixAAinPCA[,numbPCA1]+0.02, matrixAAinPCA[,numbPCA2]), 
               colour="white", alpha=1/2)+
    geom_point(size=6, aes(matrixAAinPCA[,numbPCA1], matrixAAinPCA[,numbPCA2]), 
               colour=veccol, alpha=1/2)+
    geom_text(aes(matrixAAinPCA[,numbPCA1]+0.02-1, c(matrixAAinPCA[,numbPCA2])), 
              colour=c("black"),
              label = vecAA, fontface = c("bold"), hjust=-0.25, size=5)+
    geom_text(aes(matrixAAinPCA[,numbPCA1]-0.02-1, c(matrixAAinPCA[,numbPCA2])), 
              colour=c("white"),
              label = vecAA, fontface = c("bold"), hjust=-0.25, size=5)+
    geom_text(aes(matrixAAinPCA[,numbPCA1]-1, c(matrixAAinPCA[,numbPCA2])), 
              colour=veccol,
              label = vecAA, fontface = c("bold"), hjust=-0.25, size=5)+
    labs(color="c_m quantiles clustering")+
    xlab(label = paste0("PC", numbPCA1, " ", round(eigenval[1],digits = 2) , "% var"))+
    ylab(label = paste0("PC", numbPCA2, " ", round(eigenval[2],digits = 2) , "% var"))+ 
    theme(
      panel.background = element_rect(fill = "azure3",
                                      colour = "azure3",
                                      linewidth = 0.5, linetype = "solid"),
      panel.grid.minor = element_blank())
  return(list(plotPCA))
}


treeforpPCA <- function(dataext){
  dataext$SuperClass <- rep("Fish", length(dataext$Species))
  nb         <- dim(dataext)[1]
  frm        = ~SuperClass/Class/Order/Family/Genus/Species
  phylo      <- c()
  phylo      = as.data.frame(dataext[which(names(dataext) %in% c("SuperClass", "Class", "Order", "Family", "Genus", "Species"))], stringsAsFactors = TRUE)
  phylo$Species        = stringr::str_replace(phylo$Species, " ", "_")
  phylo      <- mutate_if(phylo, is.character, as.factor)
  phylo_tree_med             = as.phylo(x = frm, data = phylo, collapse = FALSE, use.labels = TRUE)
  phylo_tree_med$edge.length = rep(1,length(phylo_tree_med$edge))

  P = phylo_tree_med
  
  if (any(is.na(P))) {
    cat("P has NA values - not allowed for")
  } else {
    cat("P looks OK") }
  return(P)
}



#############################
############################
# PLOTS FOR CROSS VALIDATION (CV)

# function to create a dataset with the cross validation outputs : data expected/data infered by phylosem/exected-infered
checkphylosemdata <- function(fileID, semID, trait, name, sample, maxCV){
  
  if (name != "all"){
    traitname = trait
    traittotest = trait
  }
  if (name == "all") {    
    traittotest = trait
    traitname = c("TOT")	
  }
  
  whichNA   <- sample[floor((1+floor(maxCV)*(fileID-1))):(floor(maxCV)*fileID)]
  if (maxCV!=floor(maxCV) & fileID==nbCV){whichNA<-c(whichNA, sample[length(sample)])}
  
  rest      <- dataext_traits[whichNA,]
  data      <- read.csv(paste0(path_phylosem_out, "/output", fileID, "_", modelname[semID], "psemFINAL", name, traitname, ".csv"))
  
  # change units of the invariants traits (Beverton and Holt)
  if(traittotest == "tm"){
    rest$tm <- log(exp(rest$tm)/rest$M)
    data$tm <- log(exp(data$tm)/data$M)
  }
  if(traittotest == "Lm"){
    rest$Lm <- rest$Lm+rest$Loo
    data$Lm <- data$Lm+data$Loo
  }
  
  #identify data to compare
  IDnonarest <- which(!is.na(rest[, traittotest]))
  species    <- rownames(rest[IDnonarest,])
  IDphylosem <- which(data$label %in% species)
  
  #extract data
  dataphylosem <- data[IDphylosem, c("label", traittotest)]
  names(dataphylosem) <- c("label", paste0(traittotest, "phylosem"))
  datafishbase <- rest[IDnonarest, c(traittotest) ]
  datafishbase <- cbind(datafishbase, rownames(rest[IDnonarest,]))
  colnames(datafishbase) <- c(traittotest, "label")
  
  dataphylosemextract <- dataphylosem[which(dataphylosem$label %in% datafishbase[,2]),]
  
  datacompare  <- full_join(dataphylosemextract, as.data.frame(datafishbase), by="label")
  datacompare[,traittotest] <- as.numeric(datacompare[,traittotest])
  
  #create comparaison column
  datacompare <- cbind(datacompare, ((datacompare[, traittotest]) - datacompare[, paste0(traittotest , "phylosem")] )*100/(datacompare[, traittotest]))
  colnames(datacompare) <- c(colnames(datacompare)[1:3], paste0(traittotest, "errorpercent"))
  
  return(datacompare)
}





plot_checkphylosemdata <- function(semID, trait, name, sample, maxCV, names_var, plot_layout=plt){
  
  for (sem in semID){
    cat("\n---------> sem", sem)
    k=0
    plotid = 0
    myplotlist <- list()
    
    for (j in trait){
      plotid = plotid+1
      names_var_i <- names_var[plotid]
      cat("\n---------> trait", j)
      if(length(sample) == length(trait)) {k=k+1}
      else {k=1}
      label    <- c()
      infered  <- c()
      expected <- c()
      errorpercent    <- c()
      
      if (checkmate::checkFileExists(paste0(path_CV,  "/dataseterror", modelname[sem], name, j, ".csv"), access = "r")==TRUE){
        dataseterrorE <- read.csv2( paste0(path_CV,  "/dataseterror", modelname[sem], name, j, ".csv"))
      }

      else{
        for (i in 1:nbCV){
          cat("\n---------> CV", i)
          toterror <- checkphylosemdata(fileID=i, semID=sem, trait=j, name=name, sample=sample[[k]], maxCV=maxCV[[k]])
          label   <- c(label, toterror[,1])
          infered <- c(infered, toterror[,2])
          expected<- c(expected, as.numeric(toterror[,3]))
          errorpercent   <- c(errorpercent, toterror[,4])
        }
        
        # http://r-statistics.co/Outlier-Treatment-With-R.html
        # identifying outliers of the dataset to mark them clearly on the plots
        dataframe_outliers <- data.frame(infered, expected)
        mod <- lm(expected ~ infered, data=dataframe_outliers)
        cooksd <- cooks.distance(mod)
        Booleand_outliers = as.numeric(cooksd>(4*mean(cooksd, na.rm=T)))
        #
        missing_cooksd <- c(0)
        
        if (length(which(is.nan(dataframe_outliers$infered)))>0){
          missing_cooksd <-which(is.nan(dataframe_outliers$infered))
          label <- label[-missing_cooksd]
          infered<-infered[-missing_cooksd]
          expected<-expected[-missing_cooksd]
          errorpercent<-errorpercent[-missing_cooksd]
        }
        
        model <- modelname[sem]
        if (sem==1){ dataseterrorE <- data.frame(label, infered, expected, errorpercent, Booleand_outliers)}
        if (sem==2){ dataseterrorM <- data.frame(label, infered, expected, errorpercent, Booleand_outliers)
        dataseterrorE <- dataseterrorM}
        if (sem==1){  write.csv2(dataseterrorE, paste0(path_CV,  "/dataseterror", modelname[sem], name, j, ".csv"))}
        if (sem==2){  write.csv2(dataseterrorM, paste0(path_CV,  "/dataseterror", modelname[sem], name, j, ".csv"))}
        
      }
      
      # Estimate percentage variance explained
      dataseterrorE <- na.omit(dataseterrorE)
      PVEe = 1- ((sum((dataseterrorE$expected-dataseterrorE$infered)^2))/
                   (sum((dataseterrorE$expected -
                           rep(times=length(dataseterrorE$expected), mean(dataseterrorE$expected)))^2)))
      ROC = roc(dataseterrorE$expected, dataseterrorE$infered)
      AUC = round(auc(ROC)*100)
      
      # plots
      
      if (sum(dataseterrorE$expected %in% c(1,0,NA))==length(dataseterrorE$expected)){
        plt = ggplot(dataseterrorE, aes(d = expected, m = infered)) + geom_roc(n.cuts = 0) + 
          ggtitle(paste0(names_var_i, "\nAUC =", AUC, "%"))+
          geom_abline (slope=1, intercept = 0, linetype = "dashed", color="Red")+
          theme_classic()+
          theme(axis.text.y = element_text(size = 9),
                axis.text.x = element_text(size = 9),
                plot.title = element_text(size=12))
      }
      if (trait == "c_m"){ 
        plt <- ggplot(dataseterrorE, aes(x = `infered`, y = `expected`)) +
          geom_point(pch = 19, alpha = 0.6, colour = ifelse(dataseterrorE$Booleand_outliers, "darkorange", "black")) +
          geom_abline(slope = 1, intercept = 0, aes(color = "red"), 
                      color = "red", linetype = "dashed", lwd = 1.5) +
          geom_smooth(aes(color = "blue")) +
          geom_text_repel(aes(label = ifelse(dataseterrorE$Booleand_outliers, label, "")),
                          size = 3.5) +
          theme_classic() +
          theme(
            plot.title = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12)) + 
          ggtitle(paste0(names_var_i, paste0("\nPVE=", label_percent()(PVEe)))) +
          scale_color_identity(name = "Model fit",
                               breaks = c("red", "blue"),
                               labels = c("x=y", "Loess"),
                               guide = "legend")
      }
      if(!(sum(dataseterrorE$expected %in% c(1,0,NA))==length(dataseterrorE$expected)) && !(trait == "c_m")) {
        dataseterrorE$LineType <- factor("LOESS")
        plt <- ggplot(dataseterrorE, aes(x=`infered`, y=`expected`)) +
          geom_point(pch=21, alpha=0.5)+ 
          geom_abline(slope = 1, intercept = 0, aes(color="red"), color="red", linetype = "dashed", lwd = 1.5)  +
          geom_smooth(aes(color = "blue", )) +
          theme_classic()+
          theme(
            plot.title = element_text(size=12),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9))+
          ggtitle(paste0(names_var_i, paste0("\nPVE=", label_percent()(PVEe))))+
          scale_color_identity(name = "Model fit",
                               breaks = c("red", "blue"),
                               labels = c("x=y", "Loess"),
                               guide = "legend") + 
          expand_limits(x = 0, y = 0)
      }
      grDevices::pdf(file=paste0(pathoutput_CV,  "/plot_CrossValidation/", modelname[sem], name, j, "plotCV.pdf"))
      print(plt)
      dev.off()
      
      plt <- plt + rremove("xlab") + rremove("ylab")
      myplotlist[[plotid]] <- plt
      
      grDevices::pdf(file=paste0(path_CV,  "/boxplot/", modelname[sem], name, j, "boxplotCV.pdf"))
      plt2 <- ggplot(dataseterrorE, aes(y=errorpercent)) + geom_boxplot()
      print(plt2)
      dev.off()
    }
    if (length(trait)>1){
      length(myplotlist)      
      par(oma=c(1, 1, 1, 1))
      par(mar=c(5, 5, 4, 2) + 0.1) 
      arranged <- ggarrange(plotlist=myplotlist, ncol=ncol, nrow=nrow,
                            labels=NULL, align="hv", 
                            font.label = list(size=7.5, color="black", face="bold", family=NULL, legend="right"), 
                            common.legend =  TRUE) +
        theme(plot.margin = margin(t = 15, r = 10, b = 15, l = 10))
      arranged <- as.list(arranged)
      if (length(arranged[[1]])==length(arranged[[2]])) {arranged_annot <-  lapply(arranged, plot_arrange)}
      else {arranged_annot <-  plot_arrange(arranged)}
      
      grDevices::pdf(file=paste0(pathoutput_CV, "/plot_CrossValidation/", modelname[sem], name, "plotarrangedCV.pdf"), height=8, width=15.37)
      print(arranged_annot)
      dev.off()
    }
  }
  return(plt)
}


plot_arrange <- function(y){
  a <- annotate_figure(y, left=text_grob("Expected trait value", rot=90, vjust=1),
                       bottom=text_grob("Infered trait value", vjust=1))
    theme(plot.margin = margin(r=15 , b = 15))
  return(a)
}


