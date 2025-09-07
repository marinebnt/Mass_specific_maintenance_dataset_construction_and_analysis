
# packages
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
library(fishtree)
library(rotl)
library(rsample)


osmosespnames <- c("Alosa alosa",                "Alosa fallax",                "Anguilla anguilla",           "Argyrosomus regius",          "Aristaeomorpha foliacea",
                   "Aristeus antennatus",         "Atherina boyeri",             "Auxis rochei",          "Belone belone",               "Boops boops",
                   "Caranx crysos",               "Chelidonichthys lucerna",     "Coris julis",                 "Coryphaena hippurus",         "Crangon crangon",
                   "Crystallogobius linearis",    "Dentex dentex",               "Dentex gibbosus",             "Dentex maroccanus",           "Dicentrarchus labrax",
                   "Diplodus annularis",          "Diplodus cervinus",           "Diplodus puntazzo",           "Diplodus sargus",       "Diplodus vulgaris",
                   "Eledone cirrhosa",            "Engraulis encrasicolus",      "Epinephelus aeneus",          "Epinephelus marginatus",      "Etrumeus sadina",
                   "Eutrigla gurnardus",          "Galeus melastomus",           "Gobius niger",                "Halobatrachus didactylus",    "Illex coindetii",
                   "Lepidorhombus whiffiagonis",  "Chelon auratus",                 "Chelon ramada",                 "Chelon saliens",                "Loligo vulgaris",
                   "Lophius budegassa",           "Lophius piscatorius",         "Merlangius merlangus",        "Merluccius merluccius",       "Micromesistius poutassou",
                   "Mugil cephalus",              "Mullus barbatus",             "Mullus surmuletus",           "Mustelus mustelus",           "Nephrops norvegicus",
                   "Octopus vulgaris",            "Pagellus acarne",             "Pagellus erythrinus",         "Pagrus pagrus",               "Palaemon serratus",
                   "Palinurus elephas",           "Parapenaeus longirostris",    "Penaeus kerathurus",          "Phycis phycis",               "Platichthys flesus",
                   "Pleuronectes platessa",       "Pomatomus saltatrix",         "Pomatoschistus marmoratus",   "Pomatoschistus minutus",      "Rhinobatos rhinobatos",
                   "Sarda sarda",                 "Sardina pilchardus",          "Sardinella aurita",           "Saurida undosquamis",         "Sciaena umbra",
                   "Scomber japonicus",          "Scomber scombrus",            "Scophthalmus maximus",        "Scorpaena notata",            "Scyliorhinus canicula",
                   "Sepia officinalis",           "Seriola dumerili",            "Serranus atricauda",          "Solea solea",                 "Sparus aurata",
                   "Sphyraena sphyraena",         "Sphyraena viridensis",        "Spicara maena",               "Spicara smaris",              "Spondyliosoma cantharus",
                   "Sprattus sprattus",           "Squilla mantis",              "Stephanolepis diaspros",      "Thunnus alalunga",            "Thunnus thynnus",
                   "Trachurus mediterraneus",     "Trachurus picturatus",        "Trachurus trachurus",         "Trachyrincus scabrus",        "Trigla lyra",
                   "Trisopterus luscus",          "Trisopterus minutus",         "Upeneus moluccensis",         "Xiphias gladius",             "Gobius ophiocephalus",
                   "Ost euphausiids")

###########################
###########################



psem <- function(semID, trait, nameCV){

  if (length(trait)==1){
    traitNAME = trait
    totvar <- var(na.omit(dataext_traits[,trait]))
    totmed <- median(na.omit(dataext_traits[,trait]))
    IDNA   <- which(!is.na(dataext_traits[,trait]))
    maxCV  <- length(IDNA)/nbCV

    repeat {
      sample  <- sample(IDNA, replace=F)  # this is to sample only among species with trait values different from NA
      CVvar  <- c()
      CVmed  <- c()
      for (i in 1:nbCV){
        sampleID <- sample[floor((1+floor(maxCV)*(i-1))):(floor(maxCV)*i)]
        CVvar[i] <- var(dataext_traits[sampleID, trait])
        CVmed[i] <- median(dataext_traits[sampleID, trait])
      }
      if (abs(100*max(abs(CVvar))/totvar - 100)<10 ){ break } #& abs(100*max(abs(CVmed))/totmed - 100)<10
    }
  }
  if (length(trait)!=1) {
    traitNAME = c("TOT")
    maxCV  <- maxtot
    sample <- sample(IDtot, replace=F)
  }



  for (semIDi in semID){
    for (i in 1:nbCV){

      whichNA <- sample[floor((1+floor(maxCV)*(i-1))):(floor(maxCV)*i)]
      data_CV <- dataext_traits
      data_CV[whichNA,trait] <- NA
      psem = phylosem(sem = as.character(modellist[[semIDi]]),
                      data = data_CV,
                      tree = P,
                      family = c(rep("binomial", 3),  rep("fixed", 18) ),
                      estimate_ou = FALSE,
                      estimate_lambda = FALSE,
                      estimate_kappa = FALSE,
                      covs = colnames(data_CV),
                      #newtonsteps = 1
      )
      matrixCSV[semIDi, i] <-  psem$opt$convergence
      assign(paste0("psemFINAL", "_rep", i, "_sem", semIDi), psem) # "_rep", i,
      df.1 = as_phylo4d(psem)
      df.2 = as(df.1, "data.frame")
      df.3 = df.2[order(df.2$label),]
      df.4 = df.3[df.3$node.type == "tip", ]
      df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]

      if (file.exists(paste0(pathoutput, "/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))){
        file.remove(paste0(pathoutput, "/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))
      }
      write.csv(df, paste0(pathoutput, "/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCV, traitNAME, ".csv"))

      matrixCSV[semIDi, i] <- psem$opt$convergence
    }

  }

  write.csv(df, paste0(pathoutput, "/", "Convergence", i, "_", "CV", nameCV, traitNAME, ".csv"))

  return(list(sample=sample, maxCV=maxCV))
}







#############################
############################
# PLOTS

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
  data      <- read.csv(paste0(pathoutput, "/output", fileID, "_", modelname[semID], "psemFINAL", name, traitname, ".csv"))

  # correct 'rest' considering the infered data is a ratio, and not the observed data (Lm and tm traits)
  # rest$Lm <- rest$Lm+rest$Loo
  # rest$tm <- log(exp(rest$tm)/rest$M)

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





plot_checkphylosemdata <- function(semID, trait, name, sample, maxCV, names_var){


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

      if (checkmate::checkFileExists(paste0(pathoutput,  "/dataseterror", modelname[sem], name, j, ".csv"), access = "r")==TRUE){
        dataseterrorE <- read.csv2( paste0(pathoutput,  "/dataseterror", modelname[sem], name, j, ".csv"))
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

      model <- modelname[sem]
      if (sem==1){ dataseterrorE <- data.frame(label, infered, expected, errorpercent, Booleand_outliers)}
      if (sem==2){ dataseterrorM <- data.frame(label, infered, expected, errorpercent, Booleand_outliers)
      dataseterrorE <- dataseterrorM}
      if (sem==1){  write.csv2(dataseterrorE, paste0(pathoutput,  "/dataseterror", modelname[sem], name, j, ".csv"))}
      if (sem==2){  write.csv2(dataseterrorM, paste0(pathoutput,  "/dataseterror", modelname[sem], name, j, ".csv"))}

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
                # legend.title = element_blank(),# legend.position = "none",
                plot.title = element_text(size=9))
      }
      else {
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
          theme(#legend.title = element_blank(),
            plot.title = element_text(size=9),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9))+
          ggtitle(paste0(names_var_i, paste0("\nPVE=", label_percent()(PVEe))))+
          # scale_color_manual(values = c("LOESS" = "blue", "x=y" = "red")) +
          # scale_linetype_manual(values = c("LOESS" = "solid", "x=y" = "dashed"))+
          scale_color_identity(name = "Model fit",
                               breaks = c("red", "blue"),
                               labels = c("x=y", "Loess"),
                               guide = "legend") +
          expand_limits(x = 0, y = 0)+
          if (trait == "c_m" && name==nameCV[3]){
            plt <- ggplot(dataseterrorE, aes(x = `infered`, y = `expected`)) +
              geom_point(pch = 19, alpha = 0.6, colour = ifelse(dataseterrorE$Booleand_outliers, "darkorange", "black")) +
              #ifelse(abs(dataseterrorE$infered - dataseterrorE$expected) > 6.5, "darkorange", "black")) +
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
      }

      grDevices::pdf(file=paste0(pathoutput,  "/plot/", modelname[sem], name, j, "plotCV.pdf"))
      print(plt)
      dev.off()

      plt <- plt + rremove("xlab") + rremove("ylab")
      myplotlist[[plotid]] <- plt

      grDevices::pdf(file=paste0(pathoutput,  "/boxplot/", modelname[sem], name, j, "boxplotCV.pdf"))
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
                            font.label = list(size=7, color="black", face="bold", family=NULL, legend="right"),
                            common.legend =  TRUE) +
        theme(plot.margin = margin(t = 15, r = 10, b = 15, l = 10))
      arranged <- as.list(arranged)
      if (length(arranged[[1]])==length(arranged[[2]])) {arranged_annot <-  lapply(arranged, plot_arrange)}
      else {arranged_annot <-  plot_arrange(arranged)}

      grDevices::pdf(file=paste0(pathoutput, "/plot/", modelname[sem], name, "plotarrangedCV.pdf"))
      print(arranged_annot)
      dev.off()
    }
  }
}


plot_arrange <- function(y){
  a <- annotate_figure(y, left=textGrob("Expected trait value", rot=90, vjust=1,
                                        gp=gpar(cex=1.1)),
                       bottom=textGrob("Infered trait value", vjust=1,
                                       gp=gpar(cex=1.1)))+
    theme(plot.margin = margin(r=15 , b = 15))
  return(a)
}



##################
##################
#  ELSE

osmosespnames <- c("Alosa alosa",                "Alosa fallax",                "Anguilla anguilla",           "Argyrosomus regius",          "Aristaeomorpha foliacea",
                   "Aristeus antennatus",         "Atherina boyeri",             "Auxis rochei",          "Belone belone",               "Boops boops",
                   "Caranx crysos",               "Chelidonichthys lucerna",     "Coris julis",                 "Coryphaena hippurus",         "Crangon crangon",
                   "Crystallogobius linearis",    "Dentex dentex",               "Dentex gibbosus",             "Dentex maroccanus",           "Dicentrarchus labrax",
                   "Diplodus annularis",          "Diplodus cervinus",           "Diplodus puntazzo",           "Diplodus sargus",       "Diplodus vulgaris",
                   "Eledone cirrhosa",            "Engraulis encrasicolus",      "Epinephelus aeneus",          "Epinephelus marginatus",      "Etrumeus sadina",
                   "Eutrigla gurnardus",          "Galeus melastomus",           "Gobius niger",                "Halobatrachus didactylus",    "Illex coindetii",
                   "Lepidorhombus whiffiagonis",  "Chelon auratus",                 "Chelon ramada",                 "Chelon saliens",                "Loligo vulgaris",
                   "Lophius budegassa",           "Lophius piscatorius",         "Merlangius merlangus",        "Merluccius merluccius",       "Micromesistius poutassou",
                   "Mugil cephalus",              "Mullus barbatus",             "Mullus surmuletus",           "Mustelus mustelus",           "Nephrops norvegicus",
                   "Octopus vulgaris",            "Pagellus acarne",             "Pagellus erythrinus",         "Pagrus pagrus",               "Palaemon serratus",
                   "Palinurus elephas",           "Parapenaeus longirostris",    "Penaeus kerathurus",          "Phycis phycis",               "Platichthys flesus",
                   "Pleuronectes platessa",       "Pomatomus saltatrix",         "Pomatoschistus marmoratus",   "Pomatoschistus minutus",      "Rhinobatos rhinobatos",
                   "Sarda sarda",                 "Sardina pilchardus",          "Sardinella aurita",           "Saurida undosquamis",         "Sciaena umbra",
                   "Scomber japonicus",          "Scomber scombrus",            "Scophthalmus maximus",        "Scorpaena notata",            "Scyliorhinus canicula",
                   "Sepia officinalis",           "Seriola dumerili",            "Serranus atricauda",          "Solea solea",                 "Sparus aurata",
                   "Sphyraena sphyraena",         "Sphyraena viridensis",        "Spicara maena",               "Spicara smaris",              "Spondyliosoma cantharus",
                   "Sprattus sprattus",           "Squilla mantis",              "Stephanolepis diaspros",      "Thunnus alalunga",            "Thunnus thynnus",
                   "Trachurus mediterraneus",     "Trachurus picturatus",        "Trachurus trachurus",         "Trachyrincus scabrus",        "Trigla lyra",
                   "Trisopterus luscus",          "Trisopterus minutus",         "Upeneus moluccensis",         "Xiphias gladius",             "Gobius ophiocephalus",
                   "Ost euphausiids")
