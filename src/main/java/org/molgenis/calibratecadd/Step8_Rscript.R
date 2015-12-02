
### new plot in gg ###
library(ggplot2)


### summary gene plots (uses output of Step 7)
#calibcaddGenes <- read.table('E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.genesumm.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
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
#ggplot of patho MAF
calibcaddGenesOrdered <- calibcaddAllGenes[order(-calibcaddAllGenes$PathoMAFThreshold),] 
p <- ggplot() +
  geom_point(data = calibcaddGenesOrdered, aes(y=PathoMAFThreshold, x=seq(1,dim(calibcaddAllGenes)[1])), alpha = 1, colour="black") +
  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_line(colour = "black"),panel.grid.minor = element_line(colour = "gray"),panel.border = element_blank(),panel.background = element_blank()) +
  labs(title="95th percentile MAF for pathogenic variants") +
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

## single gene plots
#calibcaddVariants <- read.table('E:\\Data\\clinvarcadd\\clinvar.patho.fix.snpeff.exac.withcadd.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
calibcaddVariants <- read.table('/Users/jvelde/Desktop/clinvarcadd/clinvar.patho.fix.snpeff.exac.withcadd.tsv',header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
selectGene <- "MYH7"
calibcaddVariants.selectedGene.path <- subset(calibcaddVariants, gene == selectGene & group == "PATHOGENIC")
calibcaddVariants.selectedGene.popul <- subset(calibcaddVariants, gene == selectGene & group == "POPULATION")
plot(calibcaddVariants.selectedGene.popul$pos, calibcaddVariants.selectedGene.popul$cadd, col="blue")
points(calibcaddVariants.selectedGene.path$pos, calibcaddVariants.selectedGene.path$cadd, col="red")
abline(h = mean(calibcaddVariants.selectedGene.popul$cadd), col="blue")
abline(h = mean(calibcaddVariants.selectedGene.path$cadd), col="red")

for (selectGene in unique(calibcaddVariants$gene)) {
#selectGene <- "ATM"
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
ggsave(paste("/Users/jvelde/Desktop/clinvarcadd/website/plots/",selectGene,".png", sep=""))
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

