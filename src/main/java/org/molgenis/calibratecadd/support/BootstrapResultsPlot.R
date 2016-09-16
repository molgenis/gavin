
### bootstrap analysis result processing

lightgray <- "#cccccc"; gray <- "#999999"; orange <- "#E69F00"; skyblue <- "#56B4E9"; blueishgreen <- "#009E73"
yellow <- "#F0E442"; blue <-"#0072B2"; vermillion <- "#D55E00"; reddishpurple <- "#CC79A7"

df <- read.table("/Users/joeri/github/gavin/data/other/performancebootstrap_output_usedinpaper_r0.2.r",header=TRUE)
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
                    labels=c("CADD predictive genes (737)", "CADD less predictive genes (684)", "Scarce training data genes (766)")) +
  scale_colour_manual(values=c("black", vermillion), name="GAVIN classification", breaks=c("GAVIN", "GAVINnocal"), labels=c("Gene-specific", "Genome-wide")) +
  coord_flip()

# mann-whitney-wilcoxon test and median values
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