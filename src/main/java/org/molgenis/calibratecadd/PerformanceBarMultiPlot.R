library(ggplot2)


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


df <- data.frame() 

# contents of "toolcomparison_output.r" as produced by Step9_Validation, stored here for convenience
row <- data.frame(Tool = "GAVIN", Data = "ClinVarNew", MCC = 0.7791903302275341, TN = 1537, TP = 1409, FP = 127, FN = 242, ExpBenignAsVOUS = 4, ExpPathoAsVOUS = 37, VCG = 2077, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "ClinVarNew", MCC = 0.5711285902478755, TN = 165, TP = 235, FP = 61, FN = 46, ExpBenignAsVOUS = 1442, ExpPathoAsVOUS = 1407, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "ClinVarNew", MCC = 0.71464193233354, TN = 1240, TP = 1574, FP = 425, FN = 75, ExpBenignAsVOUS = 3, ExpPathoAsVOUS = 39, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "ClinVarNew", MCC = 0.6413630478471447, TN = 865, TP = 511, FP = 170, FN = 116, ExpBenignAsVOUS = 633, ExpPathoAsVOUS = 1061, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "ClinVarNew", MCC = 0.5987171654034165, TN = 803, TP = 449, FP = 213, FN = 91, ExpBenignAsVOUS = 652, ExpPathoAsVOUS = 1148, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "ClinVarNew", MCC = 0.5174001170355323, TN = 295, TP = 457, FP = 161, FN = 77, ExpBenignAsVOUS = 1212, ExpPathoAsVOUS = 1154, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "ClinVarNew", MCC = 0.6074601623175058, TN = 1182, TP = 1345, FP = 385, FN = 237, ExpBenignAsVOUS = 101, ExpPathoAsVOUS = 106, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "ClinVarNew", MCC = 0.41337594389101545, TN = 259, TP = 420, FP = 184, FN = 95, ExpBenignAsVOUS = 1225, ExpPathoAsVOUS = 1173, VCG = 0, TotalExpertClsf = 3356); df <- rbind(df, row)
row <- data.frame(Tool = "GAVIN", Data = "MutationTaster2", MCC = 0.904582448011837, TN = 1191, TP = 136, FP = 3, FN = 23, ExpBenignAsVOUS = 0, ExpPathoAsVOUS = 2, VCG = 192, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "MutationTaster2", MCC = 0.9264078940548345, TN = 914, TP = 91, FP = 8, FN = 5, ExpBenignAsVOUS = 272, ExpPathoAsVOUS = 65, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "MutationTaster2", MCC = 0.49230045432714176, TN = 895, TP = 156, FP = 299, FN = 5, ExpBenignAsVOUS = 0, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "MutationTaster2", MCC = 0.6098037201108543, TN = 1057, TP = 141, FP = 137, FN = 20, ExpBenignAsVOUS = 0, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "MutationTaster2", MCC = 0.5831369201299343, TN = 1021, TP = 145, FP = 168, FN = 16, ExpBenignAsVOUS = 5, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "MutationTaster2", MCC = 0.6069095298545032, TN = 998, TP = 140, FP = 134, FN = 21, ExpBenignAsVOUS = 62, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "MutationTaster2", MCC = 0.5343158394950467, TN = 913, TP = 142, FP = 198, FN = 17, ExpBenignAsVOUS = 83, ExpPathoAsVOUS = 2, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "MutationTaster2", MCC = 0.7436908769908298, TN = 1002, TP = 150, FP = 79, FN = 11, ExpBenignAsVOUS = 113, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1355); df <- rbind(df, row)
row <- data.frame(Tool = "GAVIN", Data = "UMCG_Onco", MCC = 0.30344792459193126, TN = 198, TP = 23, FP = 102, FN = 3, ExpBenignAsVOUS = 1, ExpPathoAsVOUS = 0, VCG = 201, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "UMCG_Onco", MCC = 0.09998077477632912, TN = 52, TP = 1, FP = 50, FN = 0, ExpBenignAsVOUS = 199, ExpPathoAsVOUS = 25, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "UMCG_Onco", MCC = 0.20720950706112345, TN = 108, TP = 26, FP = 192, FN = 0, ExpBenignAsVOUS = 1, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "UMCG_Onco", MCC = 0.11575868150935979, TN = 220, TP = 3, FP = 64, FN = 2, ExpBenignAsVOUS = 17, ExpPathoAsVOUS = 21, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "UMCG_Onco", MCC = 0.16495721976846447, TN = 188, TP = 4, FP = 94, FN = 0, ExpBenignAsVOUS = 19, ExpPathoAsVOUS = 22, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "UMCG_Onco", MCC = 0.137535364080384, TN = 186, TP = 4, FP = 86, FN = 1, ExpBenignAsVOUS = 29, ExpPathoAsVOUS = 21, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "UMCG_Onco", MCC = 0.28682789634979156, TN = 158, TP = 26, FP = 141, FN = 0, ExpBenignAsVOUS = 2, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "UMCG_Onco", MCC = 0.08287942480796459, TN = 155, TP = 3, FP = 109, FN = 1, ExpBenignAsVOUS = 37, ExpPathoAsVOUS = 22, VCG = 0, TotalExpertClsf = 395); df <- rbind(df, row)
row <- data.frame(Tool = "GAVIN", Data = "UMCG_Various", MCC = 0.7468971473408659, TN = 1124, TP = 142, FP = 48, FN = 32, ExpBenignAsVOUS = 4, ExpPathoAsVOUS = 0, VCG = 1207, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "UMCG_Various", MCC = 0.7106809676339672, TN = 66, TP = 35, FP = 17, FN = 1, ExpBenignAsVOUS = 1093, ExpPathoAsVOUS = 138, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "UMCG_Various", MCC = 0.5613826126051713, TN = 948, TP = 165, FP = 223, FN = 9, ExpBenignAsVOUS = 5, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "UMCG_Various", MCC = 0.6073709556933514, TN = 662, TP = 53, FP = 52, FN = 11, ExpBenignAsVOUS = 462, ExpPathoAsVOUS = 110, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "UMCG_Various", MCC = 0.6005600238511881, TN = 639, TP = 52, FP = 58, FN = 8, ExpBenignAsVOUS = 479, ExpPathoAsVOUS = 114, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "UMCG_Various", MCC = 0.465004157557056, TN = 115, TP = 48, FP = 52, FN = 9, ExpBenignAsVOUS = 1009, ExpPathoAsVOUS = 117, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "UMCG_Various", MCC = 0.5715854848124708, TN = 974, TP = 159, FP = 193, FN = 15, ExpBenignAsVOUS = 9, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "UMCG_Various", MCC = 0.5268860292015254, TN = 112, TP = 53, FP = 57, FN = 3, ExpBenignAsVOUS = 1007, ExpPathoAsVOUS = 118, VCG = 0, TotalExpertClsf = 1512); df <- rbind(df, row)
row <- data.frame(Tool = "GAVIN", Data = "VariBenchTest", MCC = 0.4758697485230565, TN = 960, TP = 418, FP = 393, FN = 92, ExpBenignAsVOUS = 24, ExpPathoAsVOUS = 0, VCG = 384, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "VariBenchTest", MCC = 0.6771479204747005, TN = 774, TP = 258, FP = 113, FN = 47, ExpBenignAsVOUS = 490, ExpPathoAsVOUS = 205, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "VariBenchTest", MCC = 0.34849416435747066, TN = 631, TP = 468, FP = 746, FN = 42, ExpBenignAsVOUS = 0, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "VariBenchTest", MCC = 0.45294312509520196, TN = 1014, TP = 388, FP = 360, FN = 122, ExpBenignAsVOUS = 3, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "VariBenchTest", MCC = 0.427260893720003, TN = 884, TP = 427, FP = 490, FN = 83, ExpBenignAsVOUS = 3, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "VariBenchTest", MCC = 0.4646165568333945, TN = 924, TP = 427, FP = 435, FN = 81, ExpBenignAsVOUS = 18, ExpPathoAsVOUS = 2, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "VariBenchTest", MCC = 0.4298173770662926, TN = 738, TP = 462, FP = 562, FN = 48, ExpBenignAsVOUS = 77, ExpPathoAsVOUS = 0, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "VariBenchTest", MCC = 0.5159695593042395, TN = 1017, TP = 413, FP = 335, FN = 93, ExpBenignAsVOUS = 25, ExpPathoAsVOUS = 4, VCG = 0, TotalExpertClsf = 1887); df <- rbind(df, row)
row <- data.frame(Tool = "GAVIN", Data = "VariBenchTraining", MCC = 0.5516818511123602, TN = 7887, TP = 5337, FP = 3278, FN = 796, ExpBenignAsVOUS = 182, ExpPathoAsVOUS = 10, VCG = 3587, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "VariBenchTraining", MCC = 0.9961184506477786, TN = 6861, TP = 3792, FP = 11, FN = 8, ExpBenignAsVOUS = 4475, ExpPathoAsVOUS = 2343, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "VariBenchTraining", MCC = 0.4229251572504361, TN = 5333, TP = 5795, FP = 6014, FN = 338, ExpBenignAsVOUS = 0, ExpPathoAsVOUS = 10, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "VariBenchTraining", MCC = 0.5054718084661627, TN = 8055, TP = 4971, FP = 3187, FN = 1155, ExpBenignAsVOUS = 105, ExpPathoAsVOUS = 17, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "VariBenchTraining", MCC = 0.45413630202462435, TN = 6991, TP = 5151, FP = 4223, FN = 904, ExpBenignAsVOUS = 133, ExpPathoAsVOUS = 88, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "VariBenchTraining", MCC = 0.5016238162633618, TN = 7537, TP = 5160, FP = 3579, FN = 937, ExpBenignAsVOUS = 231, ExpPathoAsVOUS = 46, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "VariBenchTraining", MCC = 0.419323295454974, TN = 6246, TP = 5156, FP = 4432, FN = 935, ExpBenignAsVOUS = 669, ExpPathoAsVOUS = 52, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "VariBenchTraining", MCC = 0.6905403588469213, TN = 8501, TP = 5465, FP = 2132, FN = 490, ExpBenignAsVOUS = 714, ExpPathoAsVOUS = 188, VCG = 0, TotalExpertClsf = 17490); df <- rbind(df, row)

#### ADD DERIVED COLUMNS ####
df$TotalToolClsfNonVOUS <- 0
df$TotalExpertClsfNonVOUS <- 0
df$Coverage <- 0
df$Yield <- 0
for(i in 1:nrow(df)) {
  df[i,]$TotalToolClsfNonVOUS <- df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN
  df[i,]$TotalExpertClsfNonVOUS <- df[i,]$TotalToolClsfNonVOUS + df[i,]$ExpBenignAsVOUS + df[i,]$ExpPathoAsVOUS
  df[i,]$Coverage <- df[i,]$TotalToolClsfNonVOUS / df[i,]$TotalExpertClsfNonVOUS
  df[i,]$Yield <- df[i,]$MCC * df[i,]$Coverage
}

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




### bootstrap analysis result processing

df <- read.table("/Users/jvelde/github/maven/gavin/data/other/performancebootstrap_output_usedinpaper.r",header=TRUE)
df$MCCcovadj <- as.double(as.character(df$MCCcovadj))

ggplot() + geom_boxplot(data = df, aes(x = Label, fill = Calib, y = MCCcovadj)) +
  theme_classic() +
  theme(text = element_text(size = 15),  axis.title.y=element_blank(), legend.position="none") +
  ylab(expression(~MCC[covadj]~"performance")) +
# scale_fill_manual(values=c(blueishgreen, yellow, reddishpurple)) +
# scale_fill_manual(values=c("blue", "red", "purple")) +
  scale_fill_manual(values=c("gray50", "gray90", "gray70")) +
  coord_flip() +
  scale_x_discrete(limits=c("C3_GAVINnocal","C3_GAVIN","C4_GAVINnocal","C4_GAVIN", "C1_C2_GAVINnocal", "C1_C2_GAVIN"),
  labels = c("C1_C2_GAVIN" = "Gene-specific classification\nCADD predictive genes (520)\nThresholds calibrated per gene",
             "C1_C2_GAVINnocal" = "Genome-wide classification\nCADD predictive genes (520)\nThreshold 15 / MAF 0.00474",
             "C4_GAVIN" = "Gene-specific classification\nCADD less predictive genes (660)\nThresholds calibrated per gene", 
             "C4_GAVINnocal" = "Genome-wide classification\nCADD less predictive genes (660)\nThreshold 15 / MAF 0.00474", 
             "C3_GAVIN" = "Gene-specific classification\nScarce training data genes (737)\nThresholds calibrated per gene",
             "C3_GAVINnocal"="Genome-wide classification\nScarce training data genes (737)\nThreshold 15 / MAF 0.00474"
             ))


# mann-whitney-wilcoxon test and median values
calib <- subset(df, Tool == "GAVIN")
uncalib <- subset(df, Tool == "GAVINnocal")

C1_calib <- subset(calib, Calib == "C1_C2")
C1_uncalib <- subset(uncalib, Calib == "C1_C2")
median(C1_calib$MCCcovadj)
median(C1_uncalib$MCCcovadj)
wilcox.test(C1_calib$MCCcovadj, C1_uncalib$MCCcovadj)

C4_calib <- subset(calib, Calib == "C4")
C4_uncalib <- subset(uncalib, Calib == "C4")
median(C4_calib$MCCcovadj)
median(C4_uncalib$MCCcovadj)
wilcox.test(C4_calib$MCCcovadj, C4_uncalib$MCCcovadj)

C3_calib <- subset(calib, Calib == "C3")
C3_uncalib <- subset(uncalib, Calib == "C3")
median(C3_calib$MCCcovadj)
median(C3_uncalib$MCCcovadj)
wilcox.test(C3_calib$MCCcovadj, C3_uncalib$MCCcovadj)

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


