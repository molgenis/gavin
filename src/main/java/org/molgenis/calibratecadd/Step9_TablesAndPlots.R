######################################################
# Performance benchmark result processing and plotting
######################################################

version <- "r0.2"

library(ggplot2)

# Full performance results, Supplementary Table 1 in paper
source(paste("/Users/joeri/github/gavin/data/other/step9_panels_out_",version,".R",sep=""))
df$TotalToolClsfNonVOUS <- 0
df$TotalExpertClsfNonVOUS <- 0
df$Coverage <- 0
df$MCCadj <- 0
df$Sensitivity <- 0
df$Specificity <- 0
df$PPV <- 0
df$NPV <- 0
df$ACC <- 0
for(i in 1:nrow(df)) {
  df[i,]$TotalToolClsfNonVOUS <- df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN
  df[i,]$TotalExpertClsfNonVOUS <- df[i,]$TotalToolClsfNonVOUS + df[i,]$ExpBenignAsVOUS + df[i,]$ExpPathoAsVOUS
  df[i,]$Coverage <- df[i,]$TotalToolClsfNonVOUS / df[i,]$TotalExpertClsfNonVOUS
  df[i,]$MCCadj <- df[i,]$MCC * df[i,]$Coverage
  df[i,]$Sensitivity <- df[i,]$TP/(df[i,]$TP+df[i,]$FN+df[i,]$ExpPathoAsVOUS)
  df[i,]$Specificity <- df[i,]$TN/(df[i,]$TN+df[i,]$FP+df[i,]$ExpBenignAsVOUS)
  df[i,]$PPV <- df[i,]$TP/(df[i,]$TP+df[i,]$FP)
  df[i,]$NPV <- df[i,]$TN/(df[i,]$TN+df[i,]$FN)
  df[i,]$ACC <- (df[i,]$TP+df[i,]$TN)/(df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN+df[i,]$ExpPathoAsVOUS+df[i,]$ExpBenignAsVOUS)
}
#save: write.table(df, file = "~/gavinbenchmark.tsv",row.names=FALSE, sep="\t")

# Aggregated performance stats per tool, Table 3 in paper
stats <- data.frame()
for(i in unique(df$Tool)) {
  df.sub <- subset(df, Tool == i)
  medianSens <- median(df.sub$Sensitivity)
  medianSpec <- median(df.sub$Specificity)
  medianPPV <- median(df.sub$PPV)
  medianNPV <- median(df.sub$NPV)
  row <- data.frame(Tool = i, medianSens = medianSens, medianSpec = medianSpec, medianPPV = medianPPV, medianNPV = medianNPV)
  stats <- rbind(stats, row)
}
stats


# Figure 1 in paper
df.notincgd <- subset(df, Data == "NotInCGD")
ggplot() +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_point(data = df, aes(x = Sensitivity, y = Specificity, shape = Tool, colour = Tool), size=3, stroke = 3, alpha=0.75) +
  geom_text(data = df, aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, check_overlap = TRUE) +
  geom_text(data = df.notincgd , aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, colour="red",check_overlap = TRUE) +
  scale_colour_manual(values=rainbow(14)) +
  scale_shape_manual(values = c(1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)) +
  scale_x_continuous(lim=c(0.393,1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))


# Miscellaneous numbers
df.ponp2 <- subset(df, Tool == "PONP2")
# Percentage variants not classified across total benchmark variant count
(sum(df.ponp2$TotalExpertClsfNonVOUS)-sum(df.ponp2$TotalToolClsfNonVOUS)) / sum(df.ponp2$TotalExpertClsfNonVOUS)
# GAVIN per panel
df.gavin <- subset(df, Tool == "GAVIN")
df.gavin.sub <- df.gavin[,c(1,2,5,4,6,7,9,8,16,17)]
#save: write.table(df.gavin.sub, file = "~/gavinsub.tsv",row.names=FALSE, sep="\t")
# Basic calculation of pathogenic MAF threshold, CADD means, etc.
# Uses calibration results and not the benchmark output
calibcaddAllGenes <- read.table(paste("/Users/joeri/github/gavin/data/predictions/GAVIN_calibrations_",version,".tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
mean(calibcaddAllGenes$PathoMAFThreshold, na.rm = T)
mean(calibcaddAllGenes$MeanPopulationCADDScore, na.rm = T)
mean(calibcaddAllGenes$MeanPathogenicCADDScore, na.rm = T)
mean(calibcaddAllGenes$MeanDifference, na.rm = T)
sd(calibcaddAllGenes$MeanDifference, na.rm = T)
table(calibcaddAllGenes$Category)



##################################################
# Bootstrap analysis result processingand plotting
##################################################

version <- "r0.2"

library(ggplot2)

lightgray <- "#cccccc"; gray <- "#999999"; orange <- "#E69F00"; skyblue <- "#56B4E9"; blueishgreen <- "#009E73"
yellow <- "#F0E442"; blue <-"#0072B2"; vermillion <- "#D55E00"; reddishpurple <- "#CC79A7"

bootStrapResults <- paste("/Users/joeri/github/gavin/data/other/performancebootstrap_output_usedinpaper_",version,".r",sep="")
df <- read.table(bootStrapResults,header=TRUE)
df$Acc <- as.double(as.character(df$Acc))

# Figure 2 in paper
ggplot() + annotate("text", x = 1.2:6.2, y = .75, label = rep(c("Genome-wide", "Gene-specific"), 3)) + 
  geom_boxplot(data = df, aes(x = Label, y = Acc, fill = Calib, colour=Tool), size=1.1) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  ylab("Accuracy") + xlab("GAVIN classification") +
  scale_x_discrete(limits=c("C3_GAVINnocal","C3_GAVIN","C4_GAVINnocal","C4_GAVIN", "C1_C2_GAVINnocal", "C1_C2_GAVIN"),
                   labels = c("C1_C2_GAVIN" = "",
                              "C1_C2_GAVINnocal" = "",
                              "C4_GAVIN" = "", 
                              "C4_GAVINnocal" = "", 
                              "C3_GAVIN" = "",
                              "C3_GAVINnocal"="" )) +
  scale_fill_manual(values=c(blueishgreen, yellow, skyblue), 
                    name="Selected gene group",
                    breaks=c("C1_C2", "C4", "C3"),
                    labels=c("CADD predictive genes (681)", "CADD less predictive genes (732)", "Scarce training data genes (774)")) +
  scale_colour_manual(values=c("black", vermillion), name="GAVIN classification", breaks=c("GAVIN", "GAVINnocal"), labels=c("Gene-specific", "Genome-wide")) +
  coord_flip()

# mann-whitney-wilcoxon test and median accuracy values used in paper
calib <- subset(df, Tool == "GAVIN")
uncalib <- subset(df, Tool == "GAVINnocal")

C1_calib <- subset(calib, Calib == "C1_C2")
C1_uncalib <- subset(uncalib, Calib == "C1_C2")
median(C1_calib$Acc)
median(C1_uncalib$Acc)
wilcox.test(C1_calib$Acc, C1_uncalib$Acc)

C4_calib <- subset(calib, Calib == "C4")
C4_uncalib <- subset(uncalib, Calib == "C4")
median(C4_calib$Acc)
median(C4_uncalib$Acc)
wilcox.test(C4_calib$Acc, C4_uncalib$Acc)

C3_calib <- subset(calib, Calib == "C3")
C3_uncalib <- subset(uncalib, Calib == "C3")
median(C3_calib$Acc)
median(C3_uncalib$Acc)
wilcox.test(C3_calib$Acc, C3_uncalib$Acc)


