###########################
###########################
# CROSS VALIDATION

library(TMBhelper)
library(phylosem)
library(dplyr)
library(ape)
library(ggplot2)
library(ggtree)
library(phylopath)
library(stringr)
library(scales)
library(ggrepel)
require(grid)
library(ggpubr)
library(scales)
library(plotROC)
library(pROC)
library(rsample)
library(rotl)

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


psem <- function(semID, trait, nameCV){
  
  if (nameCV == "all"){
    k_CV <- rsample::vfold_cv(as.data.frame(dataset), v=nbCV, repeats = rep)
  }
  if (nameCV %in% c("c_m", "spe")){
    data_only_trait <- dataset[-which(is.na(dataset[,trait])), ]
    totvar          <- var(data_only_trait[,trait])
    if (nameCV  == "c_m"){
      data_only_trait  <- table(data_only_trait$Genus)
    }
    repeat {
      k_CV  <- rsample::vfold_cv(as.data.frame(data_only_trait), v=nbCV, repeats = rep)
      CVvar  <- c()
      CVmed  <- c()
      for (i in 1:nbCV){
        if(nameCV  == "c_m"){
          genus    <- names(data_cv_c_m)[k_CV$splits[[i]]$in_id]
          sampleID <- which(dataset[-which(is.na(dataset[,trait])), ]$Genus %in% genus)
        } else {
          sampleID <- k_CV$splits[[i]]$in_id 
        }
        CVvar[i] <- var(dataset[-which(is.na(dataset[,trait])), ][-sampleID, trait])
        CVmed[i] <- median(dataset[-which(is.na(dataset[,trait])), ][-sampleID, trait])
      }
      if (abs(100*max(abs(CVvar))/totvar - 100)<10 ){ break } #& abs(100*max(abs(CVmed))/totmed - 100)<10
    }
  }
  
  for (semIDi in semID){
    for (i in 1:nbCV){
      if (nameCV == "all"){
        whichnotNA <- k_CV$splits[[i]]$in_id
      }
      if (nameCV != "all"){
        if(nameCV  == "c_m"){
          genus      <- names(data_cv_c_m)[k_CV$splits[[i]]$in_id]
          whichnotNA <- which(dataset$Genus %in% genus)
        } else {
          whichnotNA_trait <- k_CV$splits[[i]]$in_id 
          whichnotNA <- which(dataset$SpecCode %in% data_only_trait$SpecCode[whichnotNA_trait])
        }
      }
      data_CV <- dataset_traits
      data_CV[-whichnotNA,trait] <- NA
      psem = phylosem(sem = as.character(modellist[[semIDi]]),
                      data = data_CV,
                      tree = P,
                      family = c(rep("binomial", 3),  rep("fixed", 15) ),
                      estimate_ou = FALSE,
                      estimate_lambda = FALSE,
                      estimate_kappa = FALSE,
                      covs = colnames(data_CV)
      )
      matrixCSV[semIDi, i] <-  psem$opt$convergence
      assign(paste0("psemFINAL", "_rep", i, "_sem", semIDi), psem) # "_rep", i,
      df.1 = as_phylo4d(psem)
      df.2 = as(df.1, "data.frame")
      df.3 = df.2[order(df.2$label),]
      df.4 = df.3[df.3$node.type == "tip", ]
      df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
      
      if (file.exists(paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))){
        file.remove(paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))
      }
      write.csv(df, paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))
      
      matrixCSV[semIDi, i] <- psem$opt$convergence
    }
    
  }
  
  write.csv(df, paste0(pathoutput, "/CrossValidation/", "Convergence", i, "_", "CV", nameCV, traitNAME, ".csv"))
  
  return(list(sample=k_CV))
}







#############################
############################
# PLOT

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
  
  rest      <- dataset_traits[-whichnotNA,]
  data      <- read.csv(paste0(pathoutput, "/CrossValidation/output", fileID, "_", modelname[semID], "psemFINAL", name, traitname, ".csv"))
  
  # correct 'rest' considering the infered data is a ratio, and not the observed data (Lm and tm traits)
  # rest$Lm <- rest$Lm+rest$Loo
  # rest$tm <- log(exp(rest$tm)/rest$M)
  if(traittotest == "tm"){
    rest$tm <- rest$tm - rest$M
    data$tm <- data$tm - data$M
  }
  if(traittotest == "Lm"){
    rest$Lm <- rest$Lm + rest$Loo
    data$Lm <- data$Lm + data$Loo
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
      
      if (checkmate::checkFileExists(paste0(pathoutput,  "/CrossValidation/dataseterror", modelname[sem], name, j, ".csv"), access = "r")==TRUE){
        dataseterrorE <- read.csv2( paste0(pathoutput,  "/CrossValidation/dataseterror", modelname[sem], name, j, ".csv"))
      }
      
      # else{
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
      if (sem==1){  write.csv2(dataseterrorE, paste0(pathoutput,  "/CrossValidation/dataseterror", modelname[sem], name, j, ".csv"))}
      if (sem==2){  write.csv2(dataseterrorM, paste0(pathoutput,  "/CrossValidation/dataseterror", modelname[sem], name, j, ".csv"))}
      
      # }
      
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
          #annotate("text", x = .75, y = .25, label = paste("AUC =", AUC))+
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
          geom_point(pch=21, alpha=0.5)+ #aes(col=ifelse(abs(errorpercent)>25, "big error (>25 %)", "small error (<25 %)"))) +
          geom_abline(slope = 1, intercept = 0, aes(color="red"), color="red", linetype = "dashed", lwd = 1.5)  +
          # geom_abline(slope=1, intercept = 0, linetype = "dashed", color="Red", lwd=1.5)+
          geom_smooth(aes(color = "blue", )) +
          # labs(color = "Line Type", linetype = "Line Type")+
          # geom_text(aes(label=ifelse(abs(errorpercent)>25, label," ")),
          # hjust=0.75,vjust=0.1, size=3)+#, check_overlap = TRUE)+
          theme_classic()+
          theme(
            plot.title = element_text(size=12),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9))+
          ggtitle(paste0(names_var_i, paste0("\nPVE=", label_percent()(PVEe))))+
          # scale_color_manual(values = c("LOESS" = "blue", "x=y" = "red")) +
          # scale_linetype_manual(values = c("LOESS" = "solid", "x=y" = "dashed"))+
          scale_color_identity(name = "Model fit",
                               breaks = c("red", "blue"),
                               labels = c("x=y", "Loess"),
                               guide = "legend") + 
          expand_limits(x = 0, y = 0)
##          if (trait == "c_m" && name==nameCV[3]){
#            plt <- ggplot(dataseterrorE, aes(x = `infered`, y = `expected`)) +
#              geom_point(pch = 19, alpha = 0.6, colour = ifelse(dataseterrorE$Booleand_outliers, "darkorange", "black")) +
#              #ifelse(abs(dataseterrorE$infered - dataseterrorE$expected) > 6.5, "darkorange", "black")) +
#              geom_abline(slope = 1, intercept = 0, aes(color = "red"), 
#                          color = "red", linetype = "dashed", lwd = 1.5) +
#              geom_smooth(aes(color = "blue")) +
#              geom_text_repel(aes(label = ifelse(dataseterrorE$Booleand_outliers, label, "")),
#                              size = 3.5) +
#              theme_classic() +
#              theme(
#                plot.title = element_text(size = 12),
#                axis.text.y = element_text(size = 12),
##                axis.text.x = element_text(size = 12)) + 
#              ggtitle(paste0(names_var_i, paste0("\nPVE=", label_percent()(PVEe)))) +
#              scale_color_identity(name = "Model fit",
#                                   breaks = c("red", "blue"),
#                                   labels = c("x=y", "Loess"),
#                                   guide = "legend")
#          }
      }
      
      grDevices::pdf(file=paste0(pathoutput,  "/CrossValidation/plot/", modelname[sem], name, j, "plotCV.pdf"))
      print(plt)
      dev.off()
      
      plt <- plt + rremove("xlab") + rremove("ylab")
      myplotlist[[plotid]] <- plt
      
      grDevices::pdf(file=paste0(pathoutput,  "/CrossValidation/boxplot/", modelname[sem], name, j, "boxplotCV.pdf"))
      plt2 <- ggplot(dataseterrorE, aes(y=errorpercent)) + geom_boxplot()
      print(plt2)
      dev.off()
    }
    if (length(trait)>1){
      
      # save.image(myplotlist, paste0(pathoutput, "/", modelname[sem], name, "plotarrangedCV.RData"))#save.image(plotlist, paste0(pathoutput, "/", modelname[semID], name, "plotarrangedCV.RData"))
      
      length(myplotlist)      
      par(oma=c(1, 1, 1, 1))
      par(mar=c(5, 5, 4, 2) + 0.1) 
      arranged <- ggarrange(plotlist=myplotlist, ncol=ncol, nrow=nrow,#,  widths=c(0.25,0.25),
                            labels=NULL, align="hv", 
                            font.label = list(size=7.5, color="black", face="bold", family=NULL, legend="right"), 
                            common.legend =  TRUE) +
        theme(plot.margin = margin(t = 15, r = 10, b = 15, l = 10))
      arranged <- as.list(arranged)
      if (length(arranged[[1]])==length(arranged[[2]])) {arranged_annot <-  lapply(arranged, plot_arrange)}
      else {arranged_annot <-  plot_arrange(arranged)}
      
      grDevices::pdf(file=paste0(pathoutput, "/CrossValidation/plot/", modelname[sem], name, "plotarrangedCV.pdf"))
      print(arranged_annot)
      dev.off()
    }
  }
}


#plot_arrange <- function(y){
#  a <- annotate_figure(y, left=textGrob("Expected trait value", rot=90, vjust=1,
#                                        gp=gpar(cex=1.1)),  
#                       bottom=textGrob("Infered trait value", vjust=1,
#                                       gp=gpar(cex=1.1)))+
#    theme(plot.margin = margin(r=15 , b = 15))
#  return(a)
#}

plot_arrange <- function(y){
  a <- annotate_figure(y, left=text_grob("Expected trait value", rot=90, vjust=1),
                       bottom=text_grob("Infered trait value", vjust=1))
    theme(plot.margin = margin(r=15 , b = 15))
  return(a)
}
