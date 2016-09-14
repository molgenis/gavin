
version <- "r0.3"

# Basic calculation of pathogenic MAF threshold, CADD means
calibcaddAllGenes <- read.table(paste("/Users/joeri/github/gavin/data/predictions/GAVIN_calibrations_",version,".tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
mean(calibcaddAllGenes$PathoMAFThreshold, na.rm = T)
mean(calibcaddAllGenes$MeanPopulationCADDScore, na.rm = T)
mean(calibcaddAllGenes$MeanPathogenicCADDScore, na.rm = T)
mean(calibcaddAllGenes$MeanDifference, na.rm = T)
sd(calibcaddAllGenes$MeanDifference, na.rm = T)


# CGD panel sets
source(paste("/Users/joeri/github/gavin/data/other/step9_panels_out_",version,".R",sep=""))

#### ADD DERIVED COLUMNS ####
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


#stats per tool
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


