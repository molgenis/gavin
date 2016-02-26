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

# update this with contents of "dfrows.R" as produced by Step9_Validation

row <- data.frame(Tool = "OurTool", Data = "ClinVarNew", MCC = 0.7885561054812027, TN = 1514, TP = 1202, FP = 96, FN = 227, VB = 58, VP = 259); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "ClinVarNew", MCC = 0.5711285902478755, TN = 165, TP = 235, FP = 61, FN = 46, VB = 1442, VP = 1407); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "ClinVarNew", MCC = 0.71464193233354, TN = 1240, TP = 1574, FP = 425, FN = 75, VB = 3, VP = 39); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "ClinVarNew", MCC = 0.6413630478471447, TN = 865, TP = 511, FP = 170, FN = 116, VB = 633, VP = 1061); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "ClinVarNew", MCC = 0.5987171654034165, TN = 803, TP = 449, FP = 213, FN = 91, VB = 652, VP = 1148); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "ClinVarNew", MCC = 0.5174001170355323, TN = 295, TP = 457, FP = 161, FN = 77, VB = 1212, VP = 1154); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "ClinVarNew", MCC = 0.6074601623175058, TN = 1182, TP = 1345, FP = 385, FN = 237, VB = 101, VP = 106); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "ClinVarNew", MCC = 0.41337594389101545, TN = 259, TP = 420, FP = 184, FN = 95, VB = 1225, VP = 1173); df <- rbind(df, row)
row <- data.frame(Tool = "OurTool", Data = "MutationTaster2", MCC = 0.8930554329369674, TN = 1184, TP = 106, FP = 1, FN = 23, VB = 9, VP = 32); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "MutationTaster2", MCC = 0.9264078940548345, TN = 914, TP = 91, FP = 8, FN = 5, VB = 272, VP = 65); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "MutationTaster2", MCC = 0.49230045432714176, TN = 895, TP = 156, FP = 299, FN = 5, VB = 0, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "MutationTaster2", MCC = 0.6098037201108543, TN = 1057, TP = 141, FP = 137, FN = 20, VB = 0, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "MutationTaster2", MCC = 0.5831369201299343, TN = 1021, TP = 145, FP = 168, FN = 16, VB = 5, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "MutationTaster2", MCC = 0.6069095298545032, TN = 998, TP = 140, FP = 134, FN = 21, VB = 62, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "MutationTaster2", MCC = 0.5343158394950467, TN = 913, TP = 142, FP = 198, FN = 17, VB = 83, VP = 2); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "MutationTaster2", MCC = 0.7436908769908298, TN = 1002, TP = 150, FP = 79, FN = 11, VB = 113, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "OurTool", Data = "UMCG_Onco", MCC = 0.4339042580220965, TN = 172, TP = 20, FP = 46, FN = 3, VB = 83, VP = 3); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "UMCG_Onco", MCC = 0.09998077477632912, TN = 52, TP = 1, FP = 50, FN = 0, VB = 199, VP = 25); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "UMCG_Onco", MCC = 0.20720950706112345, TN = 108, TP = 26, FP = 192, FN = 0, VB = 1, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "UMCG_Onco", MCC = 0.11575868150935979, TN = 220, TP = 3, FP = 64, FN = 2, VB = 17, VP = 21); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "UMCG_Onco", MCC = 0.16495721976846447, TN = 188, TP = 4, FP = 94, FN = 0, VB = 19, VP = 22); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "UMCG_Onco", MCC = 0.137535364080384, TN = 186, TP = 4, FP = 86, FN = 1, VB = 29, VP = 21); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "UMCG_Onco", MCC = 0.28682789634979156, TN = 158, TP = 26, FP = 141, FN = 0, VB = 2, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "UMCG_Onco", MCC = 0.08287942480796459, TN = 155, TP = 3, FP = 109, FN = 1, VB = 37, VP = 22); df <- rbind(df, row)
row <- data.frame(Tool = "OurTool", Data = "UMCG_Various", MCC = 0.8003053047735915, TN = 1073, TP = 124, FP = 22, FN = 31, VB = 81, VP = 19); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "UMCG_Various", MCC = 0.7106809676339672, TN = 66, TP = 35, FP = 17, FN = 1, VB = 1093, VP = 138); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "UMCG_Various", MCC = 0.5613826126051713, TN = 948, TP = 165, FP = 223, FN = 9, VB = 5, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "UMCG_Various", MCC = 0.6073709556933514, TN = 662, TP = 53, FP = 52, FN = 11, VB = 462, VP = 110); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "UMCG_Various", MCC = 0.6005600238511881, TN = 639, TP = 52, FP = 58, FN = 8, VB = 479, VP = 114); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "UMCG_Various", MCC = 0.465004157557056, TN = 115, TP = 48, FP = 52, FN = 9, VB = 1009, VP = 117); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "UMCG_Various", MCC = 0.5715854848124708, TN = 974, TP = 159, FP = 193, FN = 15, VB = 9, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "UMCG_Various", MCC = 0.5268860292015254, TN = 112, TP = 53, FP = 57, FN = 3, VB = 1007, VP = 118); df <- rbind(df, row)
row <- data.frame(Tool = "OurTool", Data = "VariBenchTest", MCC = 0.5733582693705775, TN = 833, TP = 303, FP = 175, FN = 85, VB = 369, VP = 122); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "VariBenchTest", MCC = 0.6771479204747005, TN = 774, TP = 258, FP = 113, FN = 47, VB = 490, VP = 205); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "VariBenchTest", MCC = 0.34849416435747066, TN = 631, TP = 468, FP = 746, FN = 42, VB = 0, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "VariBenchTest", MCC = 0.45294312509520196, TN = 1014, TP = 388, FP = 360, FN = 122, VB = 3, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "VariBenchTest", MCC = 0.427260893720003, TN = 884, TP = 427, FP = 490, FN = 83, VB = 3, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "VariBenchTest", MCC = 0.4646165568333945, TN = 924, TP = 427, FP = 435, FN = 81, VB = 18, VP = 2); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "VariBenchTest", MCC = 0.4298173770662926, TN = 738, TP = 462, FP = 562, FN = 48, VB = 77, VP = 0); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "VariBenchTest", MCC = 0.5159695593042395, TN = 1017, TP = 413, FP = 335, FN = 93, VB = 25, VP = 4); df <- rbind(df, row)
row <- data.frame(Tool = "OurTool", Data = "VariBenchTraining", MCC = 0.681176933757379, TN = 6954, TP = 4273, FP = 1330, FN = 725, VB = 3063, VP = 1145); df <- rbind(df, row)
row <- data.frame(Tool = "PONP2", Data = "VariBenchTraining", MCC = 0.9961184506477786, TN = 6861, TP = 3792, FP = 11, FN = 8, VB = 4475, VP = 2343); df <- rbind(df, row)
row <- data.frame(Tool = "CADD", Data = "VariBenchTraining", MCC = 0.4229251572504361, TN = 5333, TP = 5795, FP = 6014, FN = 338, VB = 0, VP = 10); df <- rbind(df, row)
row <- data.frame(Tool = "PROVEAN", Data = "VariBenchTraining", MCC = 0.5054718084661627, TN = 8055, TP = 4971, FP = 3187, FN = 1155, VB = 105, VP = 17); df <- rbind(df, row)
row <- data.frame(Tool = "SIFT", Data = "VariBenchTraining", MCC = 0.45413630202462435, TN = 6991, TP = 5151, FP = 4223, FN = 904, VB = 133, VP = 88); df <- rbind(df, row)
row <- data.frame(Tool = "PolyPhen2", Data = "VariBenchTraining", MCC = 0.5016238162633618, TN = 7537, TP = 5160, FP = 3579, FN = 937, VB = 231, VP = 46); df <- rbind(df, row)
row <- data.frame(Tool = "MSC", Data = "VariBenchTraining", MCC = 0.419323295454974, TN = 6246, TP = 5156, FP = 4432, FN = 935, VB = 669, VP = 52); df <- rbind(df, row)
row <- data.frame(Tool = "Condel", Data = "VariBenchTraining", MCC = 0.6905403588469213, TN = 8501, TP = 5465, FP = 2132, FN = 490, VB = 714, VP = 188); df <- rbind(df, row)


#### DERIVED COLUMNS ####
df$NrClassified <- 0
df$Total <- 0
df$Coverage <- 0
for(i in 1:nrow(df)) {
  df[i,]$NrClassified <- df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN
  df[i,]$Total <- df[i,]$NrClassified + df[i,]$VB + df[i,]$VP
  df[i,]$Coverage <- df[i,]$NrClassified / df[i,]$Total
}


#### MULTIPLOT ####
labels = c("TN", "VB", "FP", "FN", "VP", "TP");
palette <- c(vermillion, orange, skyblue, blue, lightgray, gray)
for(i in 1:nrow(df)) {
  dfrow <- df[i,]
  counts = c(dfrow$TN, dfrow$VB, dfrow$FP, dfrow$FN, dfrow$VP, dfrow$TP)
  start = c(0,cumsum(counts)); length(start) <- length(start)-1
  stop <- start+counts
  plotdata <- data.frame(start,stop,labels);
  plot <- ggplot() + geom_rect(data = plotdata, aes(xmin = start, xmax = stop, ymin = 0, ymax = 1, fill = labels)) + scale_fill_manual(values=palette) + theme(axis.text=element_blank(),rect=element_blank(),line = element_blank(), legend.position="none") +
    ggtitle(paste(dfrow$Data, "classified by",dfrow$Tool))
  assign(paste(dfrow$Data,dfrow$Tool,"Plot",sep=""), plot)
}

multiplot(
  VariBenchTrainingOurToolPlot, VariBenchTrainingCADDPlot, VariBenchTrainingMSCPlot, VariBenchTrainingPONP2Plot, VariBenchTrainingSIFTPlot, VariBenchTrainingPolyPhen2Plot, VariBenchTrainingPROVEANPlot, VariBenchTrainingCondelPlot, 
  VariBenchTestOurToolPlot, VariBenchTestCADDPlot, VariBenchTestMSCPlot, VariBenchTestPONP2Plot, VariBenchTestSIFTPlot, VariBenchTestPolyPhen2Plot, VariBenchTestPROVEANPlot, VariBenchTestCondelPlot, 
  UMCG_VariousOurToolPlot, UMCG_VariousCADDPlot, UMCG_VariousMSCPlot, UMCG_VariousPONP2Plot, UMCG_VariousSIFTPlot, UMCG_VariousPolyPhen2Plot, UMCG_VariousPROVEANPlot, UMCG_VariousCondelPlot, 
  UMCG_OncoOurToolPlot, UMCG_OncoCADDPlot, UMCG_OncoMSCPlot, UMCG_OncoPONP2Plot, UMCG_OncoSIFTPlot, UMCG_OncoPolyPhen2Plot, UMCG_OncoPROVEANPlot, UMCG_OncoCondelPlot, 
  ClinVarNewOurToolPlot, ClinVarNewCADDPlot, ClinVarNewMSCPlot, ClinVarNewPONP2Plot, ClinVarNewSIFTPlot, ClinVarNewPolyPhen2Plot, ClinVarNewPROVEANPlot, ClinVarNewCondelPlot, 
  MutationTaster2OurToolPlot, MutationTaster2CADDPlot, MutationTaster2MSCPlot, MutationTaster2PONP2Plot, MutationTaster2SIFTPlot, MutationTaster2PolyPhen2Plot, MutationTaster2PROVEANPlot, MutationTaster2CondelPlot, 
  cols=6
)


#### MCC plot #### ,fill=factor(Type)
yieldbox <- ggplot() + geom_boxplot(data = df, fill = blueishgreen, aes(x = Tool, y = MCC*Coverage)) + theme_classic() + ggtitle("Overall yield per tool across datasets (correlation x coverage)")
mccbox <- ggplot() + geom_boxplot(data = df, fill = yellow, aes(x = Tool, y = MCC)) + theme_classic() + ggtitle("Correlation per tool across datasets (Matthews Correlation Coefficient)")
coveragebox <- ggplot() + geom_boxplot(data = df, fill = reddishpurple, aes(x = Tool, y = Coverage)) + theme_classic() + ggtitle("Coverage per tool across datasets (classified/classified+unclassified)")
databox <- ggplot() + geom_boxplot(data = df, fill = gray, aes(x = Data, y = MCC*Coverage)) + theme_classic() + ggtitle("Dataset yield differences across tools")
multiplot(yieldbox, coveragebox, mccbox, databox, cols=2)

