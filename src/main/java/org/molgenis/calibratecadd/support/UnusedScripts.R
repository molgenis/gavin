library(ggplot2)
library(scales)




#### COLOURS ####
## from: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ + lightgray ##
lightgray <- "#cccccc"; gray <- "#999999"; orange <- "#E69F00"; skyblue <- "#56B4E9"; blueishgreen <- "#009E73"
yellow <- "#F0E442"; blue <-"#0072B2"; vermillion <- "#D55E00"; reddishpurple <- "#CC79A7"
cbPalette <- c(gray, orange, skyblue, blueishgreen, yellow, blue, vermillion, reddishpurple)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# original data sets
#source("/Users/joeri/github/gavin/data/other/step9_out.R")



#bit pointless but nevertheless: stats per data set (panel)
panelstats <- data.frame()
for(i in unique(df$Data)) {
  df.sub <- subset(df, Data == i)
  medianSens <- median(df.sub$Sensitivity)
  medianSpec <- median(df.sub$Specificity)
  row <- data.frame(Data = i, medianSens = medianSens, medianSpec = medianSpec)
  panelstats <- rbind(panelstats, row)
}


## nice plot
df.notincgd <- subset(df, Data == "NotInCGD")
ggplot() +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_point(data = df, aes(x = Sensitivity, y = Specificity, shape = Tool, colour = Tool), size=3, stroke = 3, alpha=0.75) +
  geom_text(data = df, aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, check_overlap = TRUE) +
  geom_text(data = df.notincgd , aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, colour="red",check_overlap = TRUE) +
  scale_colour_manual(values=rainbow(14)) +
  # scale_colour_manual(values=c(orange, "green", skyblue, "slateblue", yellow, blue, vermillion, reddishpurple, blueishgreen,"red", "blue","brown","orange", gray)) +
  scale_shape_manual(values = c(1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)) +
  scale_x_continuous(lim=c(0.393,1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))

## tiles TPR/TNR
ggplot() +
  geom_tile(data = df, aes(x = Tool, y = Data), fill = "white", colour="black") +
  theme_bw() +
  geom_tile(data = df, aes(x = Tool, y = Data, fill = Specificity, width = Sensitivity), height= 0.9) +
  geom_text(data = df, aes(x = Tool, y = Data, label = paste("Sens:",percent(Sensitivity),"\nSpec:",percent(Specificity))), colour = "black", size=2.5) +
  scale_fill_gradient(low = "aliceblue", high = "dodgerblue") +
  ylab("CGD manifestation gene panel") +
  xlab("Prediction tool")


## ppv/npv
ggplot() +
  geom_tile(data = df, aes(x = Tool, y = Data), fill = "gray", colour="black") +
  theme_bw() +
  geom_tile(data = df, aes(x = Tool, y = Data, fill = NPV, width = PPV)) +
  geom_text(data = df, aes(x = Tool, y = Data, label = paste("PPV:",percent(PPV),"\nNPV:",percent(NPV))), colour = "black", size=2.5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("CGD manifestation gene panel") +
  xlab("Prediction tool")

# acc
ggplot() +
  geom_tile(data = df, aes(x = Tool, y = Data, fill = ACC)) +
  geom_text(data = df, aes(x = Tool, y = Data, label = percent(ACC)), colour = "black", size=4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw()



#save
write.table(df, file = "~/gavinbenchmark.tsv",row.names=FALSE, sep="\t")

#plot
labels = c("TN", "VB", "FP", "FN", "VP", "TP");
palette <- c(vermillion, orange, skyblue, blue, lightgray, gray)

#### MULTIPLOT ####
for(i in 1:nrow(df)) {
  dfrow <- df[i,]
  counts = c(dfrow$TN, dfrow$ExpBenignAsVOUS, dfrow$FP, dfrow$FN, dfrow$ExpPathoAsVOUS, dfrow$TP)
  start = c(0,cumsum(counts)); length(start) <- length(start)-1
  stop <- start+counts
  plotdata <- data.frame(start,stop,labels);
  lbl <- paste("MCC[covadj] == ", round(dfrow$Yield, 2))
  plot <- ggplot() + geom_rect(data = plotdata, aes(xmin = start, xmax = stop, ymin = 0, ymax = 1, fill = labels)) + scale_fill_manual(values=palette) + theme(axis.text=element_blank(),rect=element_blank(),line = element_blank(), legend.position="none") +
    annotate("rect", xmin = sum(counts)/2-sum(counts)/4, xmax = sum(counts)/2+sum(counts)/4, ymin = .25, ymax = .75, alpha = .5, fill="white") + theme(axis.title = element_blank()) +
    annotate("text", x = sum(counts)/2, y = 0.5, label = lbl, parse=TRUE) +
    ggtitle(paste(dfrow$Data, "classified by",dfrow$Tool))
  assign(paste(dfrow$Data,dfrow$Tool,"Plot",sep=""), plot)
}

## MULTIPLOT LAYOUT 1
multiplot(
  VariBenchTrainingGAVINPlot, VariBenchTrainingCADDPlot, VariBenchTrainingMSCPlot, VariBenchTrainingPONP2Plot, VariBenchTrainingSIFTPlot, VariBenchTrainingPolyPhen2Plot, VariBenchTrainingPROVEANPlot, VariBenchTrainingCondelPlot,
  VariBenchTestGAVINPlot, VariBenchTestCADDPlot, VariBenchTestMSCPlot, VariBenchTestPONP2Plot, VariBenchTestSIFTPlot, VariBenchTestPolyPhen2Plot, VariBenchTestPROVEANPlot, VariBenchTestCondelPlot,
  UMCG_VariousGAVINPlot, UMCG_VariousCADDPlot, UMCG_VariousMSCPlot, UMCG_VariousPONP2Plot, UMCG_VariousSIFTPlot, UMCG_VariousPolyPhen2Plot, UMCG_VariousPROVEANPlot, UMCG_VariousCondelPlot,
  UMCG_OncoGAVINPlot, UMCG_OncoCADDPlot, UMCG_OncoMSCPlot, UMCG_OncoPONP2Plot, UMCG_OncoSIFTPlot, UMCG_OncoPolyPhen2Plot, UMCG_OncoPROVEANPlot, UMCG_OncoCondelPlot,
  ClinVarNewGAVINPlot, ClinVarNewCADDPlot, ClinVarNewMSCPlot, ClinVarNewPONP2Plot, ClinVarNewSIFTPlot, ClinVarNewPolyPhen2Plot, ClinVarNewPROVEANPlot, ClinVarNewCondelPlot,
  MutationTaster2GAVINPlot, MutationTaster2CADDPlot, MutationTaster2MSCPlot, MutationTaster2PONP2Plot, MutationTaster2SIFTPlot, MutationTaster2PolyPhen2Plot, MutationTaster2PROVEANPlot, MutationTaster2CondelPlot,
  cols=6
)

## MULTIPLOT LAYOUT 2
multiplot(
  UMCG_VariousGAVINPlot, UMCG_VariousCADDPlot, UMCG_VariousMSCPlot, UMCG_VariousPONP2Plot, UMCG_VariousSIFTPlot, UMCG_VariousPolyPhen2Plot, UMCG_VariousPROVEANPlot, UMCG_VariousCondelPlot,
  UMCG_OncoGAVINPlot, UMCG_OncoCADDPlot, UMCG_OncoMSCPlot, UMCG_OncoPONP2Plot, UMCG_OncoSIFTPlot, UMCG_OncoPolyPhen2Plot, UMCG_OncoPROVEANPlot, UMCG_OncoCondelPlot,
  ClinVarNewGAVINPlot, ClinVarNewCADDPlot, ClinVarNewMSCPlot, ClinVarNewPONP2Plot, ClinVarNewSIFTPlot, ClinVarNewPolyPhen2Plot, ClinVarNewPROVEANPlot, ClinVarNewCondelPlot,
  MutationTaster2GAVINPlot, MutationTaster2CADDPlot, MutationTaster2MSCPlot, MutationTaster2PONP2Plot, MutationTaster2SIFTPlot, MutationTaster2PolyPhen2Plot, MutationTaster2PROVEANPlot, MutationTaster2CondelPlot,
  VariBenchTrainingGAVINPlot, VariBenchTrainingCADDPlot, VariBenchTrainingMSCPlot, VariBenchTrainingPONP2Plot, VariBenchTrainingSIFTPlot, VariBenchTrainingPolyPhen2Plot, VariBenchTrainingPROVEANPlot, VariBenchTrainingCondelPlot,
  VariBenchTestGAVINPlot, VariBenchTestCADDPlot, VariBenchTestMSCPlot, VariBenchTestPONP2Plot, VariBenchTestSIFTPlot, VariBenchTestPolyPhen2Plot, VariBenchTestPROVEANPlot, VariBenchTestCondelPlot,
  cols=4, layout = matrix(seq(1, 48), ncol = 4, nrow = 12, byrow=T)
)


# optional: condensed per tool, sum of all data performances
for(i in unique(df$Tool)) {
  dfrow <- subset(df, Tool == i)
  counts = c(sum(dfrow$TN), sum(dfrow$ExpBenignAsVOUS), sum(dfrow$FP), sum(dfrow$FN), sum(dfrow$ExpPathoAsVOUS), sum(dfrow$TP))
  start = c(0,cumsum(counts)); length(start) <- length(start)-1
  stop <- start+counts
  plotdata <- data.frame(start,stop,labels);
  lbl <- paste("MCC[covadj] == ", round(median(dfrow$Yield), 2))
  plot <- ggplot() + geom_rect(data = plotdata, aes(xmin = start, xmax = stop, ymin = 0, ymax = 1, fill = labels)) + scale_fill_manual(values=palette) + theme(axis.text=element_blank(),rect=element_blank(),line = element_blank(), legend.position="none") +
#    annotate("rect", xmin = sum(counts)/2-sum(counts)/4, xmax = sum(counts)/2+sum(counts)/4, ymin = .25, ymax = .75, alpha = .5, fill="white") + theme(axis.title = element_blank()) +
#    annotate("text", x = sum(counts)/2, y = 0.5, label = lbl, parse=TRUE) +
#    annotate("text", x = (sum(counts)/2.25), y = 0.7, label = "median for 6 sets", size=3.5) +
    ggtitle(paste(dfrow$Tool, ", all variants", sep=""))
  assign(paste("AllData",dfrow$Tool,"Plot",sep=""), plot)
}
multiplot(
  AllDataGAVINPlot, AllDataCondelPlot, AllDataMSCPlot, AllDataCADDPlot, AllDataPolyPhen2Plot, AllDataPROVEANPlot, AllDataPONP2Plot, AllDataSIFTPlot, cols=2
)


# "Overall yield per tool across datasets (correlation x coverage)"

#### MCC plot #### ,fill=factor(Type)
yieldbox <- ggplot() + geom_boxplot(data = df, fill = blueishgreen, aes(x = Tool, y = MCC*Coverage)) + theme_classic() + theme(axis.line.x = element_line(), axis.line.y = element_line(), plot.title = element_text(size = 12), axis.text.x = element_text(angle=45, vjust=.5)) + ggtitle(expression(paste("", MCC[covadj], " per tool across datasets")))
mccbox <- ggplot() + geom_boxplot(data = df, fill = yellow , aes(x = Tool, y = MCC)) + theme_classic() + theme(axis.line.x = element_line(), axis.line.y = element_line(), plot.title = element_text(size = 12), axis.text.x = element_text(angle=45, vjust=.5)) + ggtitle("MCC per tool across datasets")
coveragebox <- ggplot() + geom_boxplot(data = df, fill = reddishpurple, aes(x = Tool, y = Coverage)) + theme_classic() + theme(axis.line.x = element_line(), axis.line.y = element_line(), plot.title = element_text(size = 12), axis.text.x = element_text(angle=45, vjust=.5)) + ggtitle("Coverage per tool across datasets")
databox <- ggplot() + geom_boxplot(data = df, fill = "gray70", aes(x = Data, y = MCC*Coverage)) + theme_classic() + theme(axis.line.x = element_line(), axis.line.y = element_line(), plot.title = element_text(size = 12), axis.text.x = element_text(angle=45, vjust=.5)) + scale_x_discrete(breaks=c("UMCG_Onco", "UMCG_Various", "VariBenchTest", "VariBenchTraining", "MutationTaster2", "ClinVarNew"), labels=c("UMCG-Onco", "UMCG-Various", "VB-Test", "VB-Training", "MT2-Test", "ClinVarNew")) + ggtitle(expression(paste("Dataset ", MCC[covadj], " distribution across tools")))
multiplot(yieldbox, coveragebox, mccbox, databox, cols=2)


#### SUBSETS ####
df.gavin <- subset(df, Tool == "GAVIN")
df.cadd <- subset(df, Tool == "CADD")
df.msc <- subset(df, Tool == "MSC")
df.ponp2 <- subset(df, Tool == "PONP2")
df.condel <- subset(df, Tool == "Condel")
df.polyphen2 <- subset(df, Tool == "PolyPhen2")
df.provean <- subset(df, Tool == "PROVEAN")
df.sift <- subset(df, Tool == "SIFT")

df.onco <- subset(df, Data == "UMCG_Onco")
df.various <- subset(df, Data == "UMCG_Various")
df.vbtest <- subset(df, Data == "VariBenchTest")
df.vbtraining <- subset(df, Data == "VariBenchTraining")
df.mt2 <- subset(df, Data == "MutationTaster2")
df.clinvar <- subset(df, Data == "ClinVarNew")

# Yield means per tool across datasets
mean(df.gavin$Yield)
mean(df.cadd$Yield)
mean(df.msc$Yield)
mean(df.ponp2$Yield)
mean(df.condel$Yield)
mean(df.polyphen2$Yield)
mean(df.provean$Yield)
mean(df.sift$Yield)

mean(df.gavin$MCC)
mean(df.cadd$MCC)
mean(df.msc$MCC)
mean(df.ponp2$MCC)
mean(df.condel$MCC)
mean(df.polyphen2$MCC)
mean(df.provean$MCC)
mean(df.sift$MCC)

median(df.gavin$MCC)
median(df.cadd$MCC)
median(df.msc$MCC)
median(df.ponp2$MCC)
median(df.condel$MCC)
median(df.polyphen2$MCC)
median(df.provean$MCC)
median(df.sift$MCC)

median(df.gavin$Yield)
median(df.cadd$Yield)
median(df.msc$Yield)
median(df.ponp2$Yield)
median(df.condel$Yield)
median(df.polyphen2$Yield)
median(df.provean$Yield)
median(df.sift$Yield)



# naive yield vs calibration coverage
nonVousClsfInCalibGenes <- df.gavin$VCG/df.gavin$TotalToolClsfNonVOUS
plot(df.gavin$Yield ~ nonVousClsfInCalibGenes)

# relative yield gain vs calibration coverage
avgYields <- c(mean(df.clinvar$Yield), mean(df.mt2$Yield), mean(df.onco$Yield), mean(df.various$Yield), mean(df.vbtest$Yield), mean(df.vbtraining$Yield))
relYield <- df.gavin$Yield-avgYields
plot(relYield ~ nonVousClsfInCalibGenes)

# LM
lmDat <- lm(relYield ~ nonVousClsfInCalibGenes)
summary(lmDat)$r.squared
lmCoeff <- coef(lmDat)

#ggplot
yieldComp <- data.frame(matrix(c(unlist(nonVousClsfInCalibGenes), unlist(relYield)), nrow=6, ncol=2, byrow=F))
colnames(yieldComp) <- c("calibcov", "rel")
plot <- ggplot() + geom_point(data = yieldComp, aes(xmin=0,xmax=1,x = calibcov, y = rel, size=5)) +
  geom_abline(intercept = lmCoeff[1], slope = lmCoeff[2]) +
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_line(colour = "darkgray"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  xlab("Percentage of variants located in calibrated genes") +
  ylab(expression(~MCC[covadj]~"improvement of GAVIN, relative to dataset average for all tools")) +
  ggtitle(expression("Yield improves when genes are calibrated ("~R^{2}~" = 0.40)"))
plot





# bonus: density plot of the same
ggplot() + geom_density(data = df, alpha = 0.5, aes(MCCcovadj, fill = Label, colour = Label)) +
  theme_classic() +
  theme(plot.title = element_text(size = 15)) +
  ggtitle("title")


# bonus: pie-chart gene calibration overview
library(ggplot2)
geneList <- c("Highly predictive (n = 257, pval < 0.01)", "Predictive (n = 263, pval < 0.05)", "Less predictive (n = 660, pval > 0.05)", "Scarce data (n = 737, pval > 0.05 with <5 samples)","'Effect impact predictive' (n = 309)", "Other (n = 829)");
df <- data.frame(
  Genes = geneList,
  value = c(257, 263, 660, 737, 309, 829)
)
df$Genes <- factor(df$Genes, levels = geneList)

bp <- ggplot(df, aes(x="", y=value, fill=Genes)) +
  geom_bar(width = 1, stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 15), axis.title.y=element_blank(), axis.title.x=element_blank()) +
  theme(legend.text=element_text(size=20), legend.title=element_text(size=25)) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank()) +
  scale_fill_manual(values=c(blue,skyblue,yellow,orange,blueishgreen,vermillion,reddishpurple)) +
  coord_polar("y", start=0)
bp







### new plot in gg ###
library(ggplot2)


### summary gene plots (uses output of Step 7)
#calibcaddGenes <- read.table('E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.genesumm.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
calibcaddVariants <- read.table('/Users/jvelde/Desktop/clinvarcadd/clinvar.patho.fix.snpeff.exac.withcadd.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
calibcaddAllGenes <- read.table('/Users/jvelde/Desktop/clinvarcadd/clinvar.patho.fix.snpeff.exac.genesumm.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
calibcaddAllGenes$PathoMAFThreshold <- as.numeric(calibcaddAllGenes$PathoMAFThreshold)
calibcaddAllGenes$MeanDifference <- as.numeric(calibcaddAllGenes$MeanDifference)
calibcaddGenes <- subset(calibcaddAllGenes, UTestPvalue != '')
#regular plot of means
plot(calibcaddGenes$MeanPopulationCADDScore, col="blue")
points(calibcaddGenes$MeanPathogenicCADDScore, col="red")
#regular density of means
plot(density(calibcaddGenes$MeanDifference))
abline(v = mean(calibcaddGenes$MeanDifference))

mean(calibcaddGenes$MeanPopulationCADDScore)
mean(calibcaddGenes$MeanPathogenicCADDScore)
mean(calibcaddGenes$MeanDifference)
sd(calibcaddGenes$MeanDifference)

#ggplot of patho MAF
calibcaddGenesOrdered <- calibcaddAllGenes[order(-calibcaddAllGenes$PathoMAFThreshold),]
p <- ggplot() +
  geom_point(data = calibcaddGenesOrdered, aes(y=PathoMAFThreshold, x=seq(1,dim(calibcaddAllGenes)[1])), alpha = 1, colour="black") +
  geom_hline(yintercept = mean(calibcaddAllGenes$PathoMAFThreshold, na.rm=T), colour = "black", size = 1, alpha = 1, linetype = 3) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_line(colour = "black"),panel.grid.minor = element_line(colour = "gray"),panel.border = element_blank(),panel.background = element_blank()) +
  labs(title=paste(sep="", "95th percentile MAFs for pathogenic variants in each gene.\nDotted line added at the mean (", format(mean(calibcaddAllGenes$PathoMAFThreshold, na.rm=T), digits=3), ")")) +
  ylab("Minor allele frequency (MAF)") +
  xlab("Genes, sorted by 95th percentile MAF value") +
  theme(legend.position = "none")
p
ggsave("~/mafplot.png", width=8, height=4.5)

##ggplot version of scatterplot
calibcaddGenesOrdered <- calibcaddGenes[order(-calibcaddGenes$MeanDifference),]
p <- ggplot() +
  geom_point(data = calibcaddGenesOrdered, aes(y=MeanPopulationCADDScore, x=seq(1,dim(calibcaddGenes)[1])), alpha = 1, colour="blue") +
  geom_point(data = calibcaddGenesOrdered, aes(y=MeanPathogenicCADDScore, x=seq(1,dim(calibcaddGenes)[1])), alpha = 1, colour="red") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title="Mean CADD scores of equalized genes: blue = population, red = pathogenic") +
  ylab("CADD score") +
  xlab("Genes by index, sorted by mean difference in pathogenic vs population mean CADD score") +
  theme(legend.position = "none")
p
ggsave("~/meansplot.png", width=8, height=4.5)

##ggplot version of density
p <- ggplot() +
  geom_line(stat="density", adjust=1, kernel = "epanechnikov", data = calibcaddGenes, size=1, alpha=1, aes(x=MeanDifference, y=..count..), colour="black") +
  geom_vline(xintercept = mean(calibcaddGenes$MeanDifference), colour = "black", size = 1, alpha = 1, linetype = "longdash") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title=paste(sep="", "Population CADD mean vs pathogenic CADD mean for equalized genes.\nMean of means = ", format(mean(calibcaddGenes$MeanDifference), digits=3))) +
  ylab("Gene density with preserved counts") +
  xlab("Gene difference of CADD scores between pathogenic mean and population mean") +
  theme(legend.position = "none")
p
ggsave("~/densitysplot.png", width=8, height=4.5)

##GG plot of all cadd scores, compared pop vs patho
calibcaddVariants.path <- subset(calibcaddVariants, group == "PATHOGENIC")
calibcaddVariants.popul <- subset(calibcaddVariants, group == "POPULATION")
p <- ggplot() +
  geom_line(stat="density", adjust=1, kernel = "epanechnikov", data = calibcaddVariants.popul, size=1, alpha=1, linetype=2, aes(x=cadd), colour="black") +
  geom_line(stat="density", adjust=1, kernel = "epanechnikov", data = calibcaddVariants.path, size=1, alpha=1, aes(x=cadd), colour="black") +
  geom_vline(colour="black", size=1, linetype=3, xintercept = 34) +
  geom_vline(colour="black", size=1, linetype=3, xintercept = 2) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_line(colour = "black"),panel.grid.minor = element_line(colour = "gray"),panel.border = element_blank(),panel.background = element_blank()) +
  labs(title="Density distribution of CADD-scored variants.\nDashed = population, solid = pathogenic. Dotted lines added at 2 and 34.") +
  ylab("Density function") +
  xlab("CADD scaled-C score") +
  theme(legend.position = "none")
p
ggsave("~/caddplot.png", width=8, height=4.5)

## single gene plots, regular R
#calibcaddVariants <- read.table('E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.withcadd.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
selectGene <- "MYH7"
calibcaddVariants.selectedGene.path <- subset(calibcaddVariants, gene == selectGene & group == "PATHOGENIC")
calibcaddVariants.selectedGene.popul <- subset(calibcaddVariants, gene == selectGene & group == "POPULATION")
plot(calibcaddVariants.selectedGene.popul$pos, calibcaddVariants.selectedGene.popul$cadd, col="blue")
points(calibcaddVariants.selectedGene.path$pos, calibcaddVariants.selectedGene.path$cadd, col="red")
abline(h = mean(calibcaddVariants.selectedGene.popul$cadd), col="blue")
abline(h = mean(calibcaddVariants.selectedGene.path$cadd), col="red")

## single gene plots, GG plot + for loop
for (selectGene in unique(calibcaddVariants$gene)) {
#selectGene <- "RYR2"
calibcaddVariants.selectedGene.path <- subset(calibcaddVariants, gene == selectGene & group == "PATHOGENIC")
calibcaddVariants.selectedGene.popul <- subset(calibcaddVariants, gene == selectGene & group == "POPULATION")
calibcaddGenes.selectedGene <- subset(calibcaddGenes, Gene == selectGene)
p <- ggplot() +
  geom_point(data = calibcaddVariants.selectedGene.popul, aes(x=pos, y=cadd), colour="blue", pch=19, alpha = .5) +
  geom_point(data = calibcaddVariants.selectedGene.path, aes(x=pos, y=cadd), colour="red", pch=19, alpha = .5) +
  geom_abline(intercept = mean(calibcaddVariants.selectedGene.popul$cadd), slope = 0, colour = "blue", size = 2, alpha = .5) +
  geom_abline(intercept = mean(calibcaddVariants.selectedGene.path$cadd), slope = 0, colour = "red", size = 2, alpha = .5) +
  geom_abline(intercept = calibcaddGenes.selectedGene$Sens95thPerCADDThreshold, slope = 0, colour = "green", size = 2, alpha = .5) +
  geom_abline(intercept = calibcaddGenes.selectedGene$Spec95thPerCADDThreshold, slope = 0, colour = "orange", size = 2, alpha = .5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.border = element_blank(),
        panel.background = element_blank()
  ) +
  labs(title=paste(selectGene, " - pathogenic mean: ", calibcaddGenes.selectedGene$MeanPathogenicCADDScore, ", population mean: ", calibcaddGenes.selectedGene$MeanPopulationCADDScore, ", p-value: ", signif(as.numeric(calibcaddGenes.selectedGene$UTestPvalue), digits=2), "\n95% sensitivity threshold (green): ", calibcaddGenes.selectedGene$Sens95thPerCADDThreshold, ", 95% specificity threshold (orange): ", calibcaddGenes.selectedGene$Spec95thPerCADDThreshold, "\nRed: ClinVar pathogenic variants - Blue: matched ExAC population variants", sep="")) +
  ylab("CADD scaled C-score") +
  xlab(paste("Genomic position [",selectGene,", chr. ",sort(unique(calibcaddVariants.selectedGene.popul$chr)),"]", sep="")) +
  theme(legend.position = "none")
#p
ggsave(paste("/Users/jvelde/Desktop/clinvarcadd/website/plots/",selectGene,".png", sep=""), width=8, height=4.5)
}

#### split variant level data for website downloads
for (selectGene in unique(calibcaddVariants$gene)) {
  #geneName <- "BRCA1"
  gendata <- subset(calibcaddVariants, gene == selectGene)
  write.table(gendata, quote=F, row.names=F, sep="\t", file=paste("/Users/jvelde/Desktop/clinvarcadd/website/data/",selectGene,".tsv", sep=""))
}



###### GG PLOT 2 - density + histogram with ggplot2, over 1 or all genes!
install.packages("ggplot2")
library(ggplot2)
for (geneName in unique(caddcalib.wilcox$gene)) {
	#geneName <- "BRCA1"
	gendata <- subset(caddcalib, gene == geneName)
	cutpoint <- caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"OptimalCutpoint"]
	palette <- c("blue","red")
	p <- ggplot() +
		geom_histogram(data = gendata, aes(x=cadd, fill=group), binwidth=1, alpha=.25, position = 'identity') +
		geom_line(stat="density", kernel = "epanechnikov", data = gendata, size=1, alpha=.75, aes(x=cadd, y=..count.., line=group, colour=group)) +
		geom_vline(colour="darkgreen", size=1, linetype="dashed", xintercept = cutpoint) +
		scale_fill_manual(values = palette) +
		scale_colour_manual(values = palette) +
		theme_bw() + theme( axis.title=element_text(size=10), plot.title=element_text(size=10), axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) +
		theme(legend.position="none") +
		xlab("CADD score for variants, blue = Matched ExAC variants (blue) vs. ClinVar pathogenic (red)") +
		ylab("Number of variants per bin (size = 1)") +
		labs(title=paste(geneName," - Youden's cutpoint: ",cutpoint,", Youden's index: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"OptimalCriterion"],2),"\nWilcoxon p-val: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"Pval"],2),", stat: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"Stat"],2),", PPV: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"PPV"],2),", NPV: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"NPV"],2),", AUC: ",signif(caddcalib.wilcox[which(caddcalib.wilcox $gene==geneName),"AUC"],2), sep=""))
	#p
	ggsave(paste("E:\\Data\\clinvarcadd\\geneplots\\",geneName,".png", sep=""))
}

##### simple, native-R one or more gene density plots
geneNames <- c("KRT10", "CHMP2B", "BUB1B", "LAMB2", "MEFV", "ANG", "PMM2", "FBN1", "RYR1", "DMD", "MYO7A", "USH2A")
for (geneName in geneNames) {
  png(filename=paste("E:\\Data\\clinvarcadd\\geneplots\\",geneName,".png",sep=""), width = 800, height = 800, pointsize = 20)
  gendata <- subset(caddcalib, gene == geneName)
  allPopul <- subset(gendata, group=="POPULATION")
  allPatho <- subset(gendata, group=="PATHOGENIC")
  c <- density(allPopul$cadd, kernel="rectangular")
  p <- density(allPatho$cadd, kernel="rectangular")
  plot(p, lwd=2, main=geneName, xlab="CADD scaled-C score (red = pathogenic, blue = common, green = cutpoint)", col="red", ylim=c(0,max(c$y, p$y)), xlim=c(0,max(allPopul$cadd,allPatho$cadd)))
  lines(c, lwd=2, col="blue")
  abline(v=caddcalib.wilcox.sorted[which(caddcalib.wilcox.sorted $gene==geneName),3], col="forestgreen", lty=3, lwd=4)
  dev.off()
}

#### Bonferoni and other facts
sign <- subset(cwcsorted, cwcsorted$Pval < 0.05/2504)
dim(sign)
notsign <- subset(cwcsorted, cwcsorted$pval >= 0.05)
dim(notsign[which(notsign$nCommon == 1 | notsign$nPatho == 1),])
notsign[which(notsign$nCommon > 5 & notsign$nPatho > 5),]

#### previous work ####
install.packages("plyr")
install.packages("OptimalCutpoints")
library(plyr)
library(OptimalCutpoints)
#### wilcox & optimal accuracy cutoff
wilcoxPlusAcc <- function(x)
{
  #wilcox
  allPopul <- subset(x, group == "POPULATION")
  allPatho <- subset(x, group == "PATHOGENIC")
  wtest <- wilcox.test(allPopul$cadd, allPatho$cadd)

  #youden (dubble neg swap?? direction = c(">"), tag.healthy = "PATHOGENIC" )
  twoclass <- subset(x, group == "PATHOGENIC" | group == "POPULATION")
  optimal.cutpoint.Youden <- optimal.cutpoints(X = "cadd", status = "group", direction = c(">"), tag.healthy = "PATHOGENIC", methods = "Youden", data = twoclass, pop.prev = NULL, categorical.cov = NULL, control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
  maxCutpointIndex <- which.max(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff)
  yoptimalcutpoint <- optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff[maxCutpointIndex]
  yppv <- unname(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$PPV[,"Value"][maxCutpointIndex])
  ynpv <- unname(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$NPV[,"Value"][maxCutpointIndex])
  yauc <- unname(optimal.cutpoint.Youden$Youden$Global$measures.acc$AUC["AUC"])
  ycriterion <- optimal.cutpoint.Youden$Youden$Global$optimal.criterion

  #return(data.frame(pval=wpval, ROC_cutoff=roc_cutoff, optimalcutpoint=max(yoptimalcutpoint), AUC=yauc, meanPatho=mean(allPatho$cadd), meanCommon=mean(allPopul$cadd), minPatho=min(allPatho$cadd), maxPatho=max(allPatho$cadd), minCommon=min(allPopul$cadd), maxCommon=max(allPopul$cadd)))
  return(data.frame(Pval=wtest$p.value, Stat=wtest$statistic, OptimalCriterion=ycriterion, OptimalCutpoint=max(yoptimalcutpoint), PPV=yppv, NPV=ynpv, AUC=yauc, nPatho=length(allPatho$cadd), nCommon=length(allPopul$cadd)))
}
caddcalib <- read.table('E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.withcadd.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
caddcalib.wilcox <- ddply(caddcalib, c("gene"), wilcoxPlusAcc)
caddcalib.wilcox.sorted = caddcalib.wilcox[ order(caddcalib.wilcox[,2]), ]
write.table(caddcalib.wilcox.sorted, row.names=F, sep="\t", file="E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.withcadd.calib.tsv", quote=F)

### HOW TO: youden index - very nice http://cran.r-project.org/web/packages/OptimalCutpoints/
gen <- subset(clin, gene == "DMD")
twoclass <- subset(gen, type == "PATHOGENIC" | type == "POPULATION")
optimal.cutpoint.Youden<-optimal.cutpoints(X = "cadd", status = "type", direction = c(">"), tag.healthy = "PATHOGENIC", methods = "Youden", data = twoclass, pop.prev = NULL, categorical.cov = NULL, control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
maxIndex <- which.max(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff)
optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff[maxIndex]
unname(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$PPV[,"Value"][maxIndex])
unname(optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$NPV[,"Value"][maxIndex])
unname(optimal.cutpoint.Youden$Youden$Global$measures.acc$AUC["AUC"])
optimal.cutpoint.Youden$Youden$Global$optimal.criterion
summary(optimal.cutpoint.Youden)

#### ROCR cutoff for one gene using accuracy... but it sucks
gen <- subset(clin, gene == "BTD")
twoclass <- subset(gen, type == "PATHOGENIC" | type == "POPULATION")
pred <- prediction(twoclass$cadd, twoclass$type)
perf <- performance(pred,measure="acc",x.measure="cutoff")
cutoff.list <- unlist(perf@x.values[[1]])
roc_cutoff <- cutoff.list[which.max(perf@y.values[[1]])]
roc_cutoff
# alternative and plot
#bestAccInd <- which.max(perf@"y.values"[[1]])
#bestMsg <- print(paste("best accuracy=",perf@"y.values"[[1]][bestAccInd], " at cutoff=",round(perf@"x.values"[[1]][bestAccInd],4),sep=""))
#plot(perf,sub=bestMsg)

###### all genes?
allPopul <- subset(datsum, type=="POPULATION")
allPatho <- subset(datsum, type=="PATHOGENIC")
#allNull <- subset(datsum, type=="NULLDIST")

c <- density(allPopul$cadd)
p <- density(allPatho$cadd)
#n <- density(allNull$cadd, na.rm=T)

plot(p, col="red", ylim=c(0,.22), xlim=c(0,50))
#lines(n, col="black")
lines(c, col="blue")

clinsum <- ddply(clin, c("gene","type"), summarize, cadd = mean(cadd))

