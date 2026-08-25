##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Multi-panel imaging analysis (PWS, widefield EU, CTCF, RNAP I/K9
# association, ChromSTEM domain radius/CVC). Large, multi-purpose
# script covering many figures across the paper:
#   - PWS metrics (Chromatin Packing Scaling/D, Chromatin Diffusion,
#     Fractional Moving Mass) across ActD/DMSO/Auxin conditions
#   - EU Signal Intensity - Widefield raincloud (Corrected Total
#     Cell Fluorescence, DMSO/6hrs ActD/1hr ActD/POLR2A 6hrs Aux)
#   - Pol 2 CTCF (Corrected Total Cell Fluorescence) degron
#     validation scatter
#   - RNAP I clustering / K9me3 association
#   - ChromSTEM Domain Radius / Domain CVC (quadrants by size category)
#
# Figure output: Fig. 1C-F (probable); Fig. 2B (confirmed --
# ggtitle("EU Signal Intensity - Widefield")); Fig. 7D-E (probable);
# Fig. 9A-C, 9E (confirmed); SI Fig. 1A-E, 1H (probable); SI Fig. 1F
# (confirmed -- ggtitle("Pol 2 CTCF")); SI Fig. 11F-G (confirmed --
# "Plot quadrants for CVC versus size", threshold lines
# geom_hline(yintercept=c(80,120))).
#
# NOTE: contains hardcoded personal macOS paths
# (/Users/lucascarter/Documents/...) alongside Quest paths -- needs
# cleanup before this is runnable by anyone else.
##################################################################

library(ggforce)
library(ggplot2)
library(PupillometryR)
library(dplyr)
library(forcats)
library(here)
library(stringr)

here()
dir.data <- here("Data/Emily/PWS") |> str_remove("Pol2 Submission/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "Pols_PWS.csv"), sep= ",", header= TRUE))
data <- data[!data$Name == "Emily_POL2RA_Untreated",] ## Scrub untreated group
data <- na.omit(data)

## Pick data to plot
dput(unique(data$Name) )
data <- data[data$Name %in% c( "Emily_POL2RA_ActD_5ugmL", "Emily_POL2RA_DMSO", "Emily_POL2RA_Auxin_6hrs", "Emily_POL1RA_Auxin_0hr", "Emily_POL1RA_Auxin_5hr"),]
data <- data[data$Name %in% c( "Emily_POL2RA_ActD_5ugmL", "Emily_POL2RA_DMSO", "Emily_POL2RA_Auxin_6hrs"),]


names <- c( "Emily_POL2RA_ActD_5ugmL", "Emily_POL2RA_DMSO", "Emily_POL2RA_Auxin_6hrs", "Emily_POL1RA_Auxin_0hr", "Emily_POL1RA_Auxin_5hr") # Get vector of experiment names
names

## Calculate P values
x <- data[data$Name ==  names[2],] 
y <- data[data$Name == names[1],]

t.test(x$D, y$D, alternative = "two.sided", var.equal = T)
t.test(x$Diffusion, y$Diffusion, alternative = "two.sided", var.equal = T)
t.test(x$FMM, y$FMM, alternative = "two.sided", var.equal = FALSE)

## Pick conditions for for-loop
data

data <- data %>% arrange(factor(Name, levels = names))

df <- data.frame()
for (i in 1:length(names)){
  
  set<- data %>%
    filter(Name == names[i],)%>% 
    mutate(factor = rep(i,time =nrow(data[data$Name == names[i],])))
  
  set <- data.frame(set$Name, set$Date, as.numeric(set$D), as.numeric(set$Diffusion), as.numeric(set$FMM), as.numeric(set$factor))
  
  df <- rbind(df,set)
  
}

colnames(df) <- c("Name","Date", "D", "Diffusion", "FMM", "factor")


## Generate separate data frames for top2a and pol data
plot <- df
plot <- na.omit(plot)
unique(plot$Name)
type(plot$factor)
##------------------------------Plotting------------------------------##

#"#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"
#"white","darkgray","#f0746e","#dc3977","#7c1d6f"

colors <- c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e")

# D Ridge plot and boxplot side by side for all Pols
plot %>%
  ggplot( aes(x=fct_reorder(Name, as.integer(factor)), y=D, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ #coord_flip() +
  #scale_fill_manual(values=c("#f0746e","#fcde9c","#525251","#dc3977","#7c1d6f"))+ # Pols
  scale_fill_manual(values=c("#009392","#9ccb86","#525251","#eeb479","#cf597e"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black","black"))+ # Pols
  #scale_fill_manual(values=c("#fcde9c","#faa476"))+ # Top2a
  #scale_color_manual(values=c("#cf597e","#cf597e","#cf597e"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  #scale_x_discrete( labels =c("TOP2A DMSO Control", "TOP2A 6 Hours Auxin"))+
  scale_x_discrete( labels =c("5 ug/mL Act D","Pol II DMSO", "Pol II 6hrs Auxin", "Pol I DMSO", "Pol I 5hrs Auxin"))+ # Pols
  scale_y_continuous(breaks = c(2.2, 2.2, 2.4, 2.6, 2.8, 3.0),limits = c(2.2, 3.0))+
  theme(legend.position="none", plot.title = element_text(size=11), text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Chromatin Packing Scaling"
          #" 40 nM Act D | 1 uM CX5461"
  ) + xlab("")

# Diffusion Ridge plot and boxplot side by side for all Pols
plot %>%
  ggplot( aes(x=fct_reorder(Name, as.integer(factor)), y=Diffusion, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ #coord_flip() +
  #scale_fill_manual(values=c("#f0746e","#fcde9c","#525251","#dc3977","#7c1d6f"))+ # Pols
  scale_fill_manual(values=c("#009392","#9ccb86","#525251","#eeb479","#cf597e"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black","black"))+ # Pols
  #scale_fill_manual(values=c("#fcde9c","#faa476"))+ # Top2a
  #scale_color_manual(values=c("#cf597e","#cf597e","#cf597e"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  #scale_x_discrete( labels =c("TOP2A DMSO Control", "TOP2A 6 Hours Auxin"))+
  scale_x_discrete( labels =c("5 ug/mL Act D","Pol II DMSO", "Pol II 6hrs Auxin", "Pol I DMSO", "Pol I 5hrs Auxin"))+ # Pols
  scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007),limits = c(0.00, 0.007))+
  theme(legend.position="none", plot.title = element_text(size=11), text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Chromatin Diffusion"
          #" 40 nM Act D | 1 uM CX5461"
  ) + xlab("")

# FMM Ridge plot and boxplot side by side for all Pols
plot %>%
  ggplot( aes(x=fct_reorder(Name, as.integer(factor)), y=FMM, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ #coord_flip() +
  #scale_fill_manual(values=c("#f0746e","#fcde9c","#525251","#dc3977","#7c1d6f"))+ # Pols
  scale_fill_manual(values=c("#009392","#9ccb86","#525251","#eeb479","#cf597e"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black","black"))+ # Pols
  #scale_fill_manual(values=c("#fcde9c","#faa476"))+ # Top2a
  #scale_color_manual(values=c("#cf597e","#cf597e","#cf597e"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  #scale_x_discrete( labels =c("TOP2A DMSO Control", "TOP2A 6 Hours Auxin"))+
  scale_x_discrete( labels =c("5 ug/mL Act D","Pol II DMSO", "Pol II 6hrs Auxin", "Pol I DMSO", "Pol I 5hrs Auxin"))+ # Pols
  scale_y_continuous(breaks = c(2e-20, 1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19, 8e-19),limits = c(0, 8e-19))+
  theme(legend.position="none", plot.title = element_text(size=11),text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Fractional Moving Mass"
          #" 40 nM Act D | 1 uM CX5461"
  ) + xlab("")

#####################################################################################
# PWS Degron analysis  Replotting of this data 042624                             
#####################################################################################

library(ggforce)
library(ggplot2)
library(PupillometryR)
library(hrbrthemes)
library(dplyr)
library(forcats)
library(here)

here()
dir.data <- here("Data/Emily/PWS") |> str_remove("Pol2 Submission/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "Pols_PWS_replot_042624.csv"), sep= ",", header= TRUE))
data <- na.omit(data)

## Pick data to plot
dput(unique(data$Name) )
names <- c("DMSO", "Emily_POL2RA_Auxin_6hrs", "Emily_POL1RA_Auxin_5hr", "ActD")# Get vector of experiment names
names <- c("DMSO", "Emily_POL2RA_Auxin_6hrs", "ActD")
names

## Calculate P values
x <- data[data$Name ==  names[2],] 
y <- data[data$Name == names[1],]

t.test(x$D, y$D, alternative = "two.sided", var.equal = T)
t.test(x$Diffusion, y$Diffusion, alternative = "two.sided", var.equal = T)
t.test(x$FMM, y$FMM, alternative = "two.sided", var.equal = FALSE)

## Pick conditions for for-loop
data

data <- data %>% arrange(factor(Name, levels = names))

df <- data.frame()
for (i in 1:length(names)){
  
  set<- data %>%
    filter(Name == names[i],)%>% 
    mutate(factor = rep(i,time =nrow(data[data$Name == names[i],])))
  
  set <- data.frame(set$Name, as.numeric(set$D), as.numeric(set$Diffusion), as.numeric(set$FMM), as.numeric(set$factor))
  
  df <- rbind(df,set)
  
}

colnames(df) <- c("Name", "D", "Diffusion", "FMM", "factor")


## Generate separate data frames for top2a and pol data
plot <- df
plot <- na.omit(plot)
unique(plot$Name)
type(plot$factor)
df <- df %>% arrange(factor(Name, levels = names))
##------------------------------Plotting------------------------------##

#"#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"
#"white","darkgray","#f0746e","#dc3977","#7c1d6f"

colors <- c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e")

# D Ridge plot and boxplot side by side for all Pols
plot %>%
  mutate(Name = fct_reorder(Name, factor)) %>% 
  ggplot( aes(x=Name, y=D, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+
  scale_fill_manual(values=c( "#009392","#9ccb86","#525251"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black"))+ # Pols
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("WT DMSO","Pol II 6hrs Auxin", "1hr Act D"))+ # Pols
  scale_y_continuous(breaks = c(2.2, 2.2, 2.4, 2.6, 2.8, 3.0),limits = c(2.2, 3.0))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16, family="Arial"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Chromatin Packing Scaling") + xlab("")

# Diffusion Ridge plot and boxplot side by side for all Pols
plot %>%
  ggplot( aes(x=fct_reorder(Name, as.integer(factor)), y=Diffusion, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ 
  scale_fill_manual(values=c("#525251", "#009392","#9ccb86","#eeb479"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black"))+ # Pols
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) + theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("WT DMSO","Pol II 6hrs Auxin", "1hr Act D"))+ # Pols
  scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003, 0.004),limits = c(0.00, 0.004))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Chromatin Diffusion") + xlab("")

# FMM Ridge plot and boxplot side by side for all Pols
plot %>%
  ggplot( aes(x=fct_reorder(Name, as.integer(factor)), y=FMM, fill=Name, colour = Name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+
  scale_fill_manual(values=c("#525251", "#009392","#9ccb86","#eeb479"))+ # Pols
  scale_color_manual(values=c("black","black","black","black","black"))+ # Pols
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +theme_classic()+
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("WT DMSO","Pol II 6hrs Auxin", "1hr Act D"))+ # Pols
  scale_y_continuous(breaks = c(2e-20, 1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19, 8e-19),limits = c(0, 8e-19))+
  theme(legend.position="none", plot.title = element_text(size=16),text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Fractional Moving Mass") + xlab("")

#####################################################################################
# PWS HCT116 ActD                               
#####################################################################################

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

dir.data <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Emily/PWS/Kants_HCT_ActD_Data"
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "010123_HCT116_ActD_D values.csv"), sep= ",", header= TRUE))
data$Name
## Calculate P values
x <- data[data$Name ==  "ActD",] 
y <- data[data$Name == "DMSO",]

t.test(x$D, y$D, alternative = "two.sided", var.equal = F)
t.test(x$Diffusion, y$Diffusion, alternative = "two.sided", var.equal = FALSE)
t.test(x$FMM, y$FMM, alternative = "two.sided", var.equal = FALSE)

#f9ddda,#f2b9c4,#e597b9,#ce78b3,#ad5fad,#834ba0,#573b88
#f3e0f7,#e4c7f1,#d1afe8,#b998dd,#9f82ce,#826dba,#63589f

# Violin plot
data %>%
  mutate(Name = fct_reorder(Name, Factor)) %>% 
  ggplot( aes(x=Name, y=D, fill=Name, colour = Name)) +
  geom_violin(trim=F, size = 0.1) +
  scale_fill_manual(values=c("#cf597e","#faa476"))+ 
  scale_color_manual(values=c("black","black","black"))+ 
  theme_classic() +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  scale_y_continuous(breaks = c(2.2, 2.4, 2.6, 2.8, 3.0),limits = c(2.2, 3.0))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Chromatin Packing Scaling"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Chromatin Packing Scaling, D") + xlab("")

data %>%
  mutate(Name = fct_reorder(Name, Factor)) %>% 
  ggplot( aes(x=Name, y=Diffusion, fill=Name, colour = Name)) +
  geom_violin(trim=F, size = 0.1) +
  scale_fill_manual(values=c("#cf597e","#faa476"))+ 
  scale_color_manual(values=c("black","black","black"))+ 
  theme_classic() +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003, 0.004),limits = c(0.00, 0.004))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Chromatin Diffusion"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Diffusion, µm^2 s^(−1)") + xlab("")

data %>%
  mutate(Name = fct_reorder(Name, Factor)) %>% 
  ggplot( aes(x=Name, y=FMM, fill=Name, colour = Name)) +
  geom_violin(trim=F, size = 0.1) +
  scale_fill_manual(values=c("#cf597e","#faa476"))+ 
  scale_color_manual(values=c("black","black","black"))+ 
  theme_classic() +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  scale_y_continuous(breaks = c(2e-20, 1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19),limits = c(0, 7e-19))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Fractional Moving Mass"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Fractional Moving Mass (g)") + xlab("")

#########################################################################################
# Pol 1 AB Widefield CTCF Results
#########################################################################################
 library(data.table)

dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Lucas/Pol 1 AB widefield 021524"
list.files(dir)

ctcf <- fread(file.path(dir, "CTCF.csv"), sep= ",", header= T)

ctcf[ctcf <0] <- 0

library(ggplot2)
library("ggExtra")

ggplot(ctcf, aes(x=ctcf$CTCF_GFP, y=ctcf$CTCF_RFP_ave, color=cond)) +theme_classic() +
  scale_color_manual(values=c("#a059a0", "#f8a07e"))+
  #scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb", "700 kb", "800 kb","900 kb", "1 MB"),breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6),limits = c(0,1e6))+ 
  #scale_x_continuous(labels = c("-1", "-0.75", "-0.5", "-0.25", "0"), breaks = c(-1, -0.75, -0.5, -0.25, 0),limits = c(-1,0.1))+ 
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position="bottom", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))+
  geom_point(size=3,alpha = 1 ) + labs(y= "Pol 1 AB", x = "Pol 1 GFP degron") 

ggplot(ctcf, aes(x=CTCF_GFP, color=cond)) + scale_color_manual(values=c("#a059a0", "#f8a07e"))+theme_classic() +
  geom_histogram(fill="white", alpha=0.5, bins = 50,position="identity") + labs(y= "Density", x = "GFP Signal") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position="bottom", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))


#########################################################################################
# Ruyi K9 STORM data
#########################################################################################

library(dplyr)
library(ggplot2)
library(PupillometryR)
library(stringr)
library(here)

here("Data/Ruyi/new_data/20240514_K9_Pol1Pol2degron (3rd)") 
# Set directory variables here
dir <- here("Data/Ruyi/new_data/20240514_K9_Pol1Pol2degron (3rd)") |> str_remove("Pol2 Submission/code")
list.files(dir)

## reorganize data
col_Data <- function(label, file, dir) {
  
  require(dplyr)  
  
  data <- read.table(file.path(dir, file), sep= ",", header= T)
  
  rownames(data) <- data$X
  data <- data[,-c(1)]
  x <- t(as.matrix(data))
  x <- data.frame(x) %>% mutate(analysis = rep(label, nrow(x)), n = seq(1,nrow(x)))
  
  return(x)
}

## gather all data REP1
#p1.int.r1 <- col_Data("pol1", "POLR1A_Auxin_interior.csv", dir=dir)
#p1.lam.r1 <- col_Data("pol1", "POLR1A_Auxin_lamin.csv", dir=dir)
#p1.nuc.r1 <- col_Data("pol1", "POLR1A_Auxin_nucleus.csv", dir=dir)
p2.int.r1 <- col_Data("pol2", "POLR2A_Auxin_interior.csv", dir=dir)
p2.lam.r1 <- col_Data("pol2", "POLR2A_Auxin_lamin.csv", dir=dir)
p2.nuc.r1 <- col_Data("pol2", "POLR2A_Auxin_nucleus.csv", dir=dir)
wt.int.r1 <- col_Data("wt", "WT_DMSO_interior.csv", dir=dir)
wt.lam.r1 <- col_Data("wt", "WT_DMSO_lamin.csv", dir=dir)
wt.nuc.r1 <- col_Data("wt", "WT_DMSO_nucleus.csv", dir=dir)

# (DMSO: N = 50; POLR1A: N=55 ; POLR2A: N= 53). 

# Set directory variables here
dir <- here("Data/Ruyi/new_data/20240522_K9_Pol1Pol2degron (4th)") |> str_remove("Pol2 Submission/code")
list.files(dir)

## gather all data REP2
#p1.int.r2 <- col_Data("pol1", "POLR1A_Auxin_interior.csv", dir=dir)
#p1.lam.r2 <- col_Data("pol1", "POLR1A_Auxin_lamin.csv", dir=dir)
#p1.nuc.r2 <- col_Data("pol1", "POLR1A_Auxin_nucleus.csv", dir=dir)
p2.int.r2 <- col_Data("pol2", "POLR2A_Auxin_interior.csv", dir=dir)
p2.lam.r2 <- col_Data("pol2", "POLR2A_Auxin_lamin.csv", dir=dir)
p2.nuc.r2 <- col_Data("pol2", "POLR2A_Auxin_nucleus.csv", dir=dir)
wt.int.r2 <- col_Data("wt", "WT_DMSO_interior.csv", dir=dir)
wt.lam.r2 <- col_Data("wt", "WT_DMSO_lamin.csv", dir=dir)
wt.nuc.r2 <- col_Data("wt", "WT_DMSO_nucleus.csv", dir=dir)

## bind
int.plt <- rbind(p1.int.r1, p2.int.r1, wt.int.r1, p1.int.r2, p2.int.r2, wt.int.r2)
lam.plt <- rbind(p1.lam.r1, p2.lam.r1, wt.lam.r1,p1.lam.r2, p2.lam.r2, wt.lam.r2)
nuc.plt <- rbind(p1.nuc.r1, p2.nuc.r1, wt.nuc.r1,p1.nuc.r2, p2.nuc.r2, wt.nuc.r2)

## bind just pol 2
int.plt <- rbind( p2.int.r1, wt.int.r1, p2.int.r2, wt.int.r2)
lam.plt <- rbind( p2.lam.r1, wt.lam.r1, p2.lam.r2, wt.lam.r2)
nuc.plt <- rbind( p2.nuc.r1, wt.nuc.r1, p2.nuc.r2, wt.nuc.r2)

colnames(int.plt ) <- c("num_blinks", "num_clust", "cluster_mean_size", "roi", "num_cluster_perUm2", "num_blink_perUm2", "analysis", "n")
colnames(lam.plt ) <- c("num_blinks", "num_clust", "cluster_mean_size", "roi", "num_cluster_perUm2", "num_blink_perUm2", "analysis", "n")
colnames(nuc.plt ) <- c("num_blinks", "num_clust", "cluster_mean_size", "roi", "num_cluster_perUm2", "num_blink_perUm2", "analysis", "n")


plot_Box <- function(data, label, ylab,col.y, col.x) {
  
  require(ggpubr)
  require(ggplot2)
  
  #my_comparisons <- list( c("pol2", "wt"), c("pol1", "wt") )
  my_comparisons <- list( c("pol2", "wt") )
  
  p <- ggplot(data, aes_string(x={{col.x}}, y = {{col.y}},  color={{col.x}}, fill={{col.x}})) +
    geom_boxplot(width = 0.65,linetype = "dashed", size=0.5, outlier.shape = 1) +
    stat_boxplot(width = 0.65,aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
    stat_boxplot(width = 0.65,geom = "errorbar", aes(ymin = ..ymax..)) + scale_color_manual(values = rep("black",3) ) +
    stat_boxplot(width = 0.65,geom = "errorbar", aes(ymax = ..ymin..)) + scale_fill_brewer(palette="Dark2") +
    ggtitle(paste0(label)) + xlab("") + ylab(ylab)+ 
    #scale_x_discrete(labels = c("Pol 1", "Pol 2", "WT"))+
    scale_x_discrete(labels = c("Pol 2", "WT"))+
    theme_classic() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0),
                            panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) + #stat_compare_means(method = "anova",label.y = 0, label.x = 0.75) +
    stat_compare_means(label = "p.signif", comparisons = my_comparisons, method = "t.test",ref.group = ".all.") 
  
  print(p)
  
}

## Blinks
plot_Box(data = int.plt, label = "Blink Density | Interior | K9 STORM", ylab = "Blinks / µm^2", col.x = "analysis", col.y="num_blink_perUm2") # interior.blinks.box
plot_Box(data = lam.plt, label = "Blink Density | Lamina | K9 STORM", ylab = "Blinks / µm^2",col.x = "analysis", col.y="num_blink_perUm2") # lamin.blinks.box
plot_Box(data = nuc.plt, label = "Blink Density | Nucleus | K9 STORM", ylab = "Blinks / µm^2",col.x = "analysis", col.y="num_blink_perUm2") # nuc.blinks.box

## Num Clusters
plot_Box(data = int.plt, label = "Cluster Density | Interior | K9 STORM", ylab = "Clusters / µm^2", col.x = "analysis", col.y="num_cluster_perUm2") # interior.clust.box
plot_Box(data = lam.plt, label = "Cluster Density | Lamina | K9 STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="num_cluster_perUm2") # lamin.clust.box
plot_Box(data = nuc.plt, label = "Cluster Density | Nucleus | K9 STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="num_cluster_perUm2") # nuc.clust.box

## cluster size
plot_Box(data = int.plt, label = "Cluster Size | Interior | K9 STORM", ylab = "Clusters / µm^2", col.x = "analysis", col.y="cluster_mean_size") # interior.clustsize.box
plot_Box(data = lam.plt, label = "Cluster Size | Lamina | K9 STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="cluster_mean_size") # lamin.clustsize.box
plot_Box(data = nuc.plt, label = "Cluster Size | Nucleus | K9 STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="cluster_mean_size") # nuc.clustsize.box

##-------------------------------------------------------## distributions of K9 cluster

library(ggplot2)
library(plyr)
library(tidyr)
library(stats)
library(forcats)

dir.r3 <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Ruyi/new_data/20240514_K9_Pol1Pol2degron (3rd)'
dir.r4 <- '/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Ruyi/new_data/20240522_K9_Pol1Pol2degron (4th)'

parse_Cluster <- function(filename, dir, cond) {
  
  require(tidyr)
  
  clus <- read.table(file.path(dir, filename), sep= ",", header= T)
  clus.out <- tidyr::gather(clus, key = cell, value = size, 2:ncol(clus))
  clus.out <- clus.out %>% dplyr::mutate(cond = paste0(cond))
  clus.out <- na.omit(clus.out)
  
  return(clus.out)
  
}

## Rep3
p1.clus <- parse_Cluster("POLR1A_Auxin_nucleus_SizeDistribution.csv", dir = dir.r3, "Pol 1")
p2.clus <- parse_Cluster("POLR2A_Auxin_nucleus_SizeDistribution.csv", dir = dir.r3, "Pol 2") 
wt.clus <- parse_Cluster("WT_DMSO_nucleus_SizeDistribution.csv", dir = dir.r3, "WT") 

## Rep4
p1.2.clus <- parse_Cluster("POLR1A_Auxin_nucleus_SizeDistribution.csv", dir = dir.r4, "Pol 1")
p2.2.clus <- parse_Cluster("POLR2A_Auxin_nucleus_SizeDistribution.csv", dir = dir.r4, "Pol 2") 
wt.2.clus <- parse_Cluster("WT_DMSO_nucleus_SizeDistribution.csv", dir = dir.r4, "WT") 

plot.df <- rbind(p1.clus, p2.clus, wt.clus, p1.2.clus, p2.2.clus, wt.2.clus)

mu <- ddply(plot.df, "cond", summarise, mean=mean(size))

## cluster density 
ggplot(plot.df, aes(x = size, fill = cond)) +
  geom_density(alpha = 0.5,color = "black")+
  scale_x_continuous(limit = c(0,5))+ scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=cond),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle(paste0("Mean: ", round(mu$mean, 2))) + xlab("RNAP I Density (1/µm^2)") + ylab("Density")

## Color scales
colors <- colorRampPalette(c("white","#CC503E","#6F4070"))(8)
colors <- c("#F0CDC7","#D36959","#A44953", "#6F4070")

library(ggridges)

## Plot distribution of all clusters
ggplot(plot.df, aes(x = size, y = cond, fill = stat(quantile))) +scale_x_continuous(limit = c(0,3))+
  stat_density_ridges(quantile_lines = F,calc_ecdf = TRUE,geom = "density_ridges_gradient") +theme_classic() +
  scale_fill_manual(values = colors)+ theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("Cluster Size") + xlab("Clusters per µm^2") + ylab("")

#########################################################################################
# Ruyi EU STORM 052024
#########################################################################################
rm(list = ls())
gc()

library(here)
library(stringr)

# Set directory variables here
dir <- here("Data/Ruyi/new_data/20240524_eu_STORM/processed_data") |> str_remove("Pol2 Submission/code")
list.files(dir)

## gather all data REP1
p2.r1.nuc <- col_Data("pol2", "POLR2A_Auxin_nucleus.csv", dir=dir) ## old pol 2 data (sparse)
p2.r2.nuc <- col_Data("pol2", "POLR2A_Auxin_nucleus_071724.csv", dir=dir) ## updated pol 2 data from 071724
wt.r1.nuc <- col_Data("wt", "WT_DMSO_nucleus.csv", dir=dir)
wt.r2.nuc <- col_Data("wt", "WT_DMSO_nucleus_071724.csv", dir=dir)
# (POLR1A: N = 22; POLR2A: N = 18; DMSO: N= 28).  
nuc.plt <- rbind(wt.r1.nuc, wt.r2.nuc, p2.r1.nuc, p2.r2.nuc)

colnames(nuc.plt ) <- c("num_blinks", "num_clust", "cluster_mean_size", "roi", "num_cluster_perUm2", "num_blink_perUm2", "analysis", "n")

plot_Box <- function(data, label, ylab,col.y, col.x) {
  
  require(ggpubr)
  require(ggplot2)
  
  my_comparisons <- list( c( "wt","pol2") )
  colors <- c("#39b185","#e88471")
  
  p <- ggplot(data, aes_string(x={{col.x}}, y = {{col.y}},  color={{col.x}}, fill={{col.x}})) +
    geom_boxplot(width = 0.75,linetype = "dashed", size=0.5, outlier.shape = 1) +
    stat_boxplot(width = 0.75,aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
    stat_boxplot(width = 0.75,geom = "errorbar", aes(ymin = ..ymax..)) +scale_color_manual(values = rep("black",3) ) +
    stat_boxplot(width = 0.75,geom = "errorbar", aes(ymax = ..ymin..)) + scale_fill_manual(values = colors) +
    ggtitle(paste0(label)) + xlab("") + ylab(ylab)+ scale_x_discrete(labels = c("Pol 2", "WT"))+
    #scale_y_continuous(breaks = c(500,1000,1500,2000), labels = c("0.5k", "1k", "1.5k", "2k"))+ # for blink density
    theme_classic() + theme(legend.position="",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=18, family="Arial"), axis.text.x = element_text(angle = 0),
                            panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) + #stat_compare_means(method = "anova",label.y = 0, label.x = 0.75) +
    stat_compare_means(label = "p.format", comparisons = my_comparisons, method = "t.test",ref.group = ".all.") 
  
  print(p)
  
}

## Blinks
plot_Box(data = nuc.plt, label = "Blink Density | Nucleus | EU STORM", ylab = "Blinks / µm^2",col.x = "analysis", col.y="num_blink_perUm2") # nuc.blinks.box
plot_Box(data = nuc.plt, label = "Cluster Density | Nucleus | EU STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="num_cluster_perUm2") # nuc.clust.box
plot_Box(data = nuc.plt, label = "Cluster Size | Nucleus | EU STORM", ylab = "Clusters / µm^2",col.x = "analysis", col.y="cluster_mean_size") # nuc.clustsize.box

## t.test comparisons
x <- p1.nuc$num_cluster_nucleus.per.um.2
y <- wt.nuc$num_cluster_nucleus.per.um.2

x <- p1.nuc$num_blink_nucleus.per.um.2
y <- wt.nuc$num_blink_nucleus.per.um.2

x <- p1.nuc$cluster_mean_size_nucleus
y <- wt.nuc$cluster_mean_size_nucleus

x <- p1.nuc$num_blink_nucleus.per.um.2
y <- wt.nuc$num_blink_nucleus.per.um.2

t.test(x, y)

#########################################################################################
# calculate probability on STORM pol 2 ActD data
#########################################################################################
library(ggplot2)
library(plyr)
library(stats)

dir <- here("Data/Nico/data_3color_STORM") |> str_remove("Pol2 Submission/code")
list.files(dir)

pol2.prob <- read.table(file.path(dir, "pol2densitycsv.csv"), sep= ",", header= T)
K9.clus <- read.table(file.path(dir, "k9clustersize.csv"), sep= ",", header= T)

pol2.prob <- na.omit(pol2.prob )
K9.clus <- na.omit(K9.clus)

## https://aakinshin.net/posts/kde-bw/ for more info on BW
## KDE plots https://www.r-bloggers.com/2019/08/arguments-of-statsdensity/
# bw = c("nrd0", "nrd", "ucv", "bcv", "SJ")
ko.p2 <- data.frame(y = pol2.prob$Act.D/sum(pol2.prob$Act.D), x = pol2.prob$Act.D, group = rep("Act D", length(pol2.prob$Act.D)))
ctrl.p2 <- data.frame(y = pol2.prob$DMSO/sum(pol2.prob$DMSO), x = pol2.prob$DMSO, group = rep("DMSO", length(pol2.prob$DMSO)))
p2.df <- rbind(ko.p2,ctrl.p2)

ko.p2 <- data.frame(y = pol2.prob$Act.D, x = pol2.prob$Act.D, group = rep("Act D", length(pol2.prob$Act.D)))
ctrl.p2 <- data.frame(y = pol2.prob$DMSO, x = pol2.prob$DMSO, group = rep("DMSO", length(pol2.prob$DMSO)))
p2.df <- rbind(ko.p2,ctrl.p2)

ko.k9 <- data.frame(y = K9.clus$Act.D/sum(K9.clus$Act.D), x = K9.clus$Act.D, group = rep("Act D", length(K9.clus$Act.D)))
ctrl.k9 <- data.frame(y = K9.clus$DMSO/sum(K9.clus$DMSO), x = K9.clus$DMSO, group = rep("DMSO", length(K9.clus$DMSO)))
k9.df <- rbind(ko.k9,ctrl.k9)

cols <- c("#cf597e","#faa476")
cols <- c("#C6011F")

mu <- ddply(p2.df, "group", summarise, mean=median(x))
mu <- ddply(ctrl.p2, "group", summarise, mean=median(x)) ## just dmso
mu
## pol2.density.plot
ggplot(ctrl.p2, aes(x = x, fill = group)) +
  geom_density(alpha = 0.8, color = NA)+
  scale_x_continuous(limit = c(0,1000))+
  scale_fill_manual(values = cols) + theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=group),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  #ggtitle("1 Hour Act D (5 µg/mL)") + 
  ggtitle("") + 
  xlab("RNAP II Density (1/µm^2)") + ylab("Density")

mu <- ddply(ctrl.k9, "group", summarise, mean=median(x))
mu
mu$mean[1]-mu$mean[2]
## pol2.density.plot
ggplot(ctrl.k9, aes(x = x, fill = group)) +
  geom_density(alpha = 0.8, color = NA) + scale_x_continuous(limit = c(0,2))+
  scale_fill_manual(values = cols) + theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=group),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  #ggtitle("1 Hour Act D (5 µg/mL)") +
  ggtitle("") +
  xlab("Relative K9me3 Cluster Area") + ylab("Density")

## Plot association data
assoc <- read.table(file.path(dir, "pol2_associationPlot.csv"), sep= ",", header= T)
assoc <- assoc[c(1,3),] # remove act d
assoc$fact <- as.numeric(c(1,2,1,2)) ## act d included
assoc$fact <- as.numeric(c(1,2)) ## just dmso
cols <- c("#cf597e","#faa476")
cols <- c("black")

assoc %>% 
  mutate(cond = fct_reorder(cond, as.numeric(fact))) %>% 
  ggplot(aes(x = loc, y = value, fill = cond)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + scale_fill_manual(values = cols) + theme_classic()+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  #ggtitle("1 Hour Act D (5 µg/mL)") + 
  ggtitle("") + 
  xlab("") + ylab("RNAPII Association") +scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label= paste0(round(value*100, 2),"%")), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5) + theme(legend.position = "none")

dir.data <- "/SMLM/H3K9Me3 Pol2 Assoc SMLM"
file.path(root,dir.data, "pol2densitycsv.csv")
pol2.prob <- read.table(file.path(root,dir.data, "pol2densitycsv.csv"), sep= ",", header= T)
#########################################################################################
# calculate probability on STORM pol 1 data
#########################################################################################
library(ggplot2)
library(plyr)
library(stats)
library(forcats)

dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Nico/K9_Pol1_Data_CTRL"
list.files(dir)

## Nico gave me these backwards for some reason. So they are entered in backwards here, to fix that
k9.clus <- read.table(file.path(dir, "Pol1_dataCompiled_pol1DensV3.csv"), sep= ",", header= T)
pol1.prob <- read.table(file.path(dir, "Pol1_dataCompiled_k9clustersizeV3.csv"), sep= ",", header= T)

pol1.prob <- na.omit(pol1.prob )
k9.clus <- na.omit(k9.clus)

## https://aakinshin.net/posts/kde-bw/ for more info on BW
## KDE plots https://www.r-bloggers.com/2019/08/arguments-of-statsdensity/
# bw = c("nrd0", "nrd", "ucv", "bcv", "SJ")


p1.df <- data.frame(x = sort(pol1.prob$CTRL), group = rep("Pol 1", nrow(pol1.prob)))
k9.df <- data.frame(x = sort(k9.clus$K9_Pol1), group = rep("K9", nrow(k9.clus) ) )

mu <- ddply(p1.df, "group", summarise, mean=mean(x))

## pol2.density.plot
ggplot(p1.df, aes(x = x, fill = "black")) +
  geom_density(alpha = 0.8, fill = "black",color = "black")+
  scale_x_continuous(limit = c(0,500))+
  scale_fill_manual(values = cols) + theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=group),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle(paste0("Mean: ", round(mu$mean, 2))) + xlab("RNAP I Density (1/µm^2)") + ylab("Density")

mu <- ddply(k9.df, "group", summarise, mean=mean(x))

## pol2.density.plot
ggplot(k9.df, aes(x = x, fill = group)) +
  geom_density(alpha = 0.8, trim = F, fill = "black",color = "black") + scale_x_continuous(limit = c(0,2))+
  scale_fill_manual(values = cols) + theme_classic()+
  geom_vline(data=mu, aes(xintercept=mean, color=group),linetype="dashed")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle(paste0("Mean: ", round(mu$mean, 2))) + xlab("Relative K9me3 Cluster Area") + ylab("Density")

## Plot association data
assoc <- read.table(file.path(dir, "Pol1_rawData_pol1AffilFreqV2.csv"), sep= ",", header= T)
colnames(assoc) <- c("X", "Unassociated", "Inside One Cluster", "Inside Mult. Clusters")
assoc <- assoc %>% gather(key=loc, value = value, 2:4) %>% mutate(cell = as.numeric(X)) %>%  dplyr::rename(factor=X) ## prepare for plotting

cols <- c("#A5AB77","#39b185","#009392")

assoc<- ddply(assoc , "cell",transform, label_ypos=cumsum(value))

assoc %>% 
  mutate(cell = fct_reorder(as.character(cell), as.numeric(factor))) %>% 
  ggplot(aes(x = cell, y = value, fill = loc)) + 
  geom_bar(stat="identity", position="stack", color = NA) + scale_fill_manual(values = cols) + theme_classic()+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("") + xlab("Cell") + ylab("RNAPI K9 Association") +scale_y_continuous(labels = scales::percent) +
  geom_text(aes(y=label_ypos, label= paste0(round(value*100, 2),"%")), vjust=1.6, color="white", size=3.5)


#########################################################################################
# Act D Confocal Coefficient of Variation Analysis
#########################################################################################

library(ggplot2)
library(beeswarm)
library(ggbeeswarm)

# Set directory variables here
dir.data <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Emily/ActD_CoV_Analysis"
list.files(dir.data)

## Read in data
df <- data.frame(read.table(file.path(dir.data, "CVanalysis_result.csv"), sep= ",", header= TRUE))

## Calculate P values
x <- data[data$Cond ==  "Untreated",] 
y <- data[data$Cond == "1hr ActD",]

t.test(x, y, alternative = "two.sided", var.equal = F)

beeswarm(df$CoV ~ df$name, pch = 19,col =c("#D8BFD8", "#5302a3") )
legend("topright", legend = c("1 Hr Act D", "DMSO"),
       col = c("#D8BFD8", "#5302a3"), pch = 19)

## Boxplot 
p <- df %>%
  mutate(name= fct_reorder(name, as.numeric(factor))) %>% 
  ggplot(aes(x= name, y=CoV, color = name)) + theme_classic() +
  geom_beeswarm(cex = 3)+
  #geom_boxplot(color = "black", notch=F, notchwidth = 0) + 
  scale_y_continuous(limits = c(0.1,  0.5), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_x_discrete( labels =c( "DMSO", "1 hr Act D"))+
  scale_color_manual(values=c("#cf597e","#faa476")) +
  labs(title=" 1 Hour Actinomycin D") +
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 0, size = 16),
        axis.text.y = element_text(angle = 0,size = 16))+
  ylab("Coefficient of Variation")
p

#########################################################################################
# Pol1 and Pol2 Widefield CTCF results
#########################################################################################
library(data.table)
library(ggplot2)
library("ggside")

dir<- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Lucas/051524_pol1_pol2_AB"

ctcf <- fread(file.path(dir, "CTCF_values.csv"), sep= ",", header= T)
ctcf[ctcf <0] <- 0 ## remove negative CTCF values

plot <- ctcf[ctcf$cond == "pol1_aux" | ctcf$cond == "pol1_dmso",]
plot <- ctcf[ctcf$cond == "pol2_aux" | ctcf$cond == "pol2_dmso",]

p<- ggplot(plot, aes(x=gfp, y=rfp, color=cond)) +theme_classic() +
  scale_color_manual(values=c("#7570B3", "#E6AB02"))+ ggtitle("Pol 2 CTCF")+
  scale_y_continuous(labels = seq(0,max(plot$gfp),1000),breaks = seq(0,max(plot$gfp),1000),limits = c(0,5500))+ 
  scale_x_continuous(labels = seq(0,max(plot$rfp),1000),breaks = seq(0,max(plot$rfp),1000),limits = c(0,max(plot$rfp)-300))+  
  geom_point(size=2,alpha = 0.8 ) + theme(legend.position="bottom") + labs(y= "AB IF CTCF", x = "Degron IF CTCF") +
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 90))
p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.5, xfill = cond), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.5, yfill = cond), position = "dodge")+
  scale_xfill_manual(values = c("#7570B3", "#E6AB02")) +
  scale_yfill_manual(values = c("#7570B3", "#E6AB02")) + theme(ggside.panel.scale = .3)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 

ggplot(ctcf, aes(x=CTCF_GFP, color=cond)) + scale_color_manual(values=c("#a059a0", "#f8a07e"))+theme_classic() +
  geom_histogram(fill="white", alpha=0.5, bins = 50,position="identity")

#########################################################################################
# CTCF values Pol1 AB and Pol2 AB degron confirmation
#########################################################################################
library(data.table)
library(ggplot2)
library("ggside")

dir<- "/Users/lucascarter/Documents/IBiS/Backman_Lab/Transcription Publication/Data/Lucas/051524_pol1_pol2_AB"

ctcf <- fread(file.path(dir, "CTCF_values.csv"), sep= ",", header= T)
ctcf[ctcf <0] <- 0 ## remove negative CTCF values

plot <- ctcf[ctcf$cond == "pol1_aux" | ctcf$cond == "pol1_dmso",]
plot <- ctcf[ctcf$cond == "pol2_aux" | ctcf$cond == "pol2_dmso",]

p<- ggplot(plot, aes(x=gfp, y=rfp, color=cond)) +theme_classic() +
  scale_color_manual(values=c("#7570B3", "#E6AB02"))+ ggtitle("Pol 2 CTCF")+
  scale_y_continuous(labels = seq(0,max(plot$gfp),1000),breaks = seq(0,max(plot$gfp),1000),limits = c(0,5500))+ 
  scale_x_continuous(labels = seq(0,max(plot$rfp),1000),breaks = seq(0,max(plot$rfp),1000),limits = c(0,max(plot$rfp)-300))+  
  geom_point(size=2,alpha = 0.8 ) + theme(legend.position="bottom") + labs(y= "AB IF CTCF", x = "Degron IF CTCF") +
  theme( plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 90))
p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.5, xfill = cond), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.5, yfill = cond), position = "dodge")+
  scale_xfill_manual(values = c("#7570B3", "#E6AB02")) +
  scale_yfill_manual(values = c("#7570B3", "#E6AB02")) + theme(ggside.panel.scale = .3)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 

ggplot(ctcf, aes(x=CTCF_GFP, color=cond)) + scale_color_manual(values=c("#a059a0", "#f8a07e"))+theme_classic() +
  geom_histogram(fill="white", alpha=0.5, bins = 50,position="identity")

#########################################################################################
# EU Click It Widefield CTCF results
#########################################################################################

library(ggplot2)
library(tidyverse)
library(PupillometryR)
library(here)
library(stringr)
here()

# Set directory variables here
dir.data <- here("Data/Lucas/EU widefield 110923/CTCF_data") |> str_remove("Pol2 Submission/code")
list.files(dir.data)

## Read in data
df <- data.frame(read.table(file.path(dir.data, "CTCF_data_final.csv"), sep= ",", header= TRUE))

## Calculate P values
x <- df[df$condition ==  "Pol2_6hrs_Aux_R1",]
x <- df[df$condition ==  "wt_6hractd",]
x <- df[df$condition ==  "wt_1hractd",]
# (POLR1A: N = 172; POLR2A: N = 142; 6hrActD: N=20; 1hrActD: N=28; DMSO: N= 204)
y <- df[df$condition == "all_dmso",]

t.test(as.numeric(x$CTCF), as.numeric(y$CTCF), alternative = "two.sided", var.equal = F)
unique(df$condition)

colors <- c("#9ccb86","#e9e29c","#eeb479","#e88471")
colors <- c("#009392","#39b185","#9ccb86","#e9e29c")

unique(df$condition)
df <- df[!df$condition == "rad21_6hrs_Aux_R1" & !df$condition =="Pol1_6hrs_Aux_R1",] ## Remove Rad21 ## Is actually old Rad21 with bad auxin

df %>%
  mutate(condition= fct_reorder(condition, as.numeric(factor))) %>% 
  ggplot( aes(x=condition, y=as.numeric(CTCF), fill=condition, colour = condition)) +
  geom_flat_violin(color = NA,trim = F, position = position_nudge(x = .15))+ 
  scale_fill_manual(values=colors)+ # Top2a
  scale_color_manual(values=c("black","black","black", "black","black","black","black"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey"), notch = FALSE) +
  theme_classic() +
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c( "Pol 2 6hrs Aux","6hrs Act D", "1hr Act D", "DMSO Control"))+
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000),limits = c(-5000, 50000))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 0)) + coord_flip() +
  ggtitle("EU Signal Intensity - Widefield"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Corrected Total Cell Fluorescence") +xlab("") +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
                                                                 axis.text.y = element_text(angle = 0,size = 12))



#####################################################################################
# ChromSTEM Act D re-plotting                              
#####################################################################################

library(ggforce)
library(PupillometryR)

dir.data <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/bench_projects/Lab Members Results/Wing Shun"
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "Chrom_STEM_Domains.csv"), sep= ",", header= TRUE))

## Calculate P values
x <- data[data$Cond ==  "ActD",] 
y <- data[data$Cond == "DMSO",]

t.test(x$D, y$D, alternative = "two.sided", var.equal = F)
t.test(x$Rf, y$Rf, alternative = "two.sided", var.equal = FALSE)
t.test(x$CVC, y$CVC, alternative = "two.sided", var.equal = FALSE)

# D Ridge plot and boxplot side by side for all Pols
data %>%
  mutate(Cond = fct_reorder(Cond, factor)) %>% 
  ggplot( aes(x=Cond, y=D, fill=Cond, colour = Cond)) +
  geom_flat_violin(color = NA, trim = F, position = position_nudge(x = .15))+
  scale_fill_manual(values=c("#cf597e","#faa476"))+ # Top2a
  scale_color_manual(values=c("black","black","black"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic() +
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  scale_y_continuous(breaks = c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0),limits = c(2.0, 3.0))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Chromatin Packing Scaling"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Chromatin Packing Scaling, D")

# domain radius
data %>%
  mutate(Cond = fct_reorder(Cond, factor)) %>% 
  ggplot( aes(x=Cond, y=Rf, fill=Cond, colour = Cond)) +
  geom_flat_violin(color = NA,trim = F, position = position_nudge(x = .15))+ 
  scale_fill_manual(values=c("#cf597e","#faa476"))+ # Top2a
  scale_color_manual(values=c("black","black","black"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic() +
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  #scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007),limits = c(0.00, 0.007))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Domain Radius"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Domain Radius (nm)")

# CVC
data %>%
  mutate(Cond = fct_reorder(Cond, factor)) %>% 
  ggplot( aes(x=Cond, y=CVC, fill=Cond, colour = Cond)) +
  geom_flat_violin(color = NA,trim = F, position = position_nudge(x = .15))+ 
  scale_fill_manual(values=c("#cf597e","#faa476"))+ # Top2a
  scale_color_manual(values=c("black","black","black"))+ # Top2a
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic() +
  geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+
  scale_x_discrete( labels =c("DMSO Control", "1 Hour Act D"))+
  #scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007),limits = c(0.00, 0.007))+
  theme(legend.position="none", plot.title = element_text(size=16), text = element_text(size = 16), axis.text.x = element_text(angle = 0)) +
  ggtitle("Domain CVC"
          #" 40 nM Act D | 1 uM CX5461"
  ) + ylab("Chromatin Volume Concentration")

#####################################################################################
# ChromSTEM Pol1 plotting                              
#####################################################################################

library(ggforce)
library(PupillometryR)
library(here)
library(tidyverse)

# Set directory variables here
dir.data <- here("Data/Wing Shun/Pol1 Chrom STEM") |> str_remove("Pol2 Submission/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "alldomains_081725.csv"), sep= ",", header= TRUE))

brewer.pal(6, "Dark2")
cols <- c("#7570B3", "#E6AB02")

## Function to plot chromSTEM data
plot_violin <- function(data, title, xlab, ylab, y.pos,col.y, col.x, cols, comps) {
  
  require(ggpubr)
  require(ggplot2)
  
  my_comparisons <- comps
  
  p <- ggplot(data, aes_string(x={{col.x}}, y = {{col.y}},  color={{col.x}}, fill={{col.x}})) +
    geom_flat_violin(color = NA,trim = F, position = position_nudge(x = .15))+ 
    stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+ scale_color_manual(values=rep("black", length(unique(data$cond))+1))+ scale_fill_manual(values = cols) +
    ggtitle(paste0(title)) + xlab("") + ylab(ylab)+ scale_x_discrete(labels = xlab)+
    theme_classic() + theme(legend.position="",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0),
                            panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) + 
    stat_compare_means(#label = "p.format",
      label = "p.signif",
      comparisons = my_comparisons, 
      #method = "t.test", 
      ref.group = ".all.", ## for multiple comparisons
      label.y = y.pos
    ) 
  
  print(p)
  
}
# labels: "WT_DMSO"       "Pol1_6hrs_Aux" "Pol1_DMSO"     "Pol2_6hrs_Aux" "Pol2_DMSO" 

data.p <- data[data$cond =="Pol1_6hrs_Aux" | data$cond =="WT_DMSO" ,]
data.p <- data[data$cond =="Pol1_6hrs_Aux" | data$cond =="Pol1_DMSO" ,]
data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="WT_DMSO" ,]
data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="Pol2_DMSO" ,]

cond1 = "Pol2_6hrs_Aux"
cond2 = "Pol2_DMSO"

data.p <- data.p[data.p$D >1,]
## single comparisons
plot_violin(data=data.p, title = "Chromatin Packing Scaling", ylab = "D", xlab=c("Pol 2", "WT"), y.pos = 2.95, col.y = "D", col.x="cond", cols= cols, comps = list(c(cond1, cond2 )))
plot_violin(data=data.p, title = "Domain Radius", ylab = "Domain Radius (nm)", xlab=c("Pol 2", "WT"), y.pos = 225, col.y = "Rmax", col.x="cond", cols= cols, list(c(cond1, cond2  )))
plot_violin(data=data.p, title = "Domain CVC", ylab = "Chromatin Volume Concentration", xlab=c("Pol 2", "WT"), y.pos = 0.5, col.y = "Density", col.x="cond", cols= cols, comps = list(c(cond1, cond2  ))) 

cols <- c("#7570B3", "#E6AB02","#66A61E")

## Multiple comparisons
plot_violin(data=data, title = "Chromatin Packing Scaling", ylab = "D", xlab=c("Pol 1 Aux","Pol 1 DMSO",  "WT DMSO"), y.pos = c(2.9, 2.95, 3), col.y = "D", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" ))) 
plot_violin(data=data, title = "Domain Radius", ylab = "Domain Radius (nm)", xlab=c("Pol 1 Aux","Pol 1 DMSO", "WT DMSO"), y.pos = c(210, 233, 245), col.y = "Rmax", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" ))) 
plot_violin(data=data, title = "Domain CVC", ylab = "Chromatin Volume Concentration", xlab=c("Pol 1 Aux","Pol 1 DMSO", "WT DMSO"), y.pos = c(0.45,0.5,0.55), col.y = "Density", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" )))

##-------------------------------------------------------## Plot quadrants


library(ggside)
library(ggpubr)

cond1 = "Pol2_6hrs_Aux"
cond2 = "Pol2_DMSO"

data.p <- data[data$cond ==cond1 | data$cond ==cond2 ,]
data.p$lab <- NA
pe <- median(data.p$PackEff[data.p$cond ==cond2])
rad <- median(data.p$Rmax[data.p$cond ==cond2])

paste0("packing efficiency for ", cond2, ": ", mean(data.p$PackEff[data.p$cond ==cond2]) )
paste0("domain size for ", cond2, ": ",mean(data.p$Rmax[data.p$cond ==cond2]) )
paste0("scaling for ", cond2, ": ",mean(data.p$D[data.p$cond ==cond2]) )

paste0("packing efficiency for ", cond1, ": ", mean(data.p$PackEff[data.p$cond ==cond1]) )
paste0("domain size for ", cond1, ": ",mean(data.p$Rmax[data.p$cond ==cond1]) )
paste0("scaling for ", cond1, ": ",mean(data.p$D[data.p$cond ==cond1]) )
##-- reorder data if needed for plotting layers

# Reorder rows by 'Name' column alphabetically from Z to A
data.p.sorted <- data.p[order(data.p$cond, decreasing = TRUE), ]
data.p <- data.p.sorted 

##-- end reorder

## Split data up into quarter
data.p$lab[data.p$PackEff <= pe & data.p$Rmax >= rad] <- "Decaying"
data.p$lab[data.p$PackEff >= pe & data.p$Rmax >= rad] <- "Large Mature"
data.p$lab[data.p$PackEff >= pe & data.p$Rmax <= rad] <- "Small Mature"
data.p$lab[data.p$PackEff <= pe & data.p$Rmax <= rad] <- "Nascent"
#data.p$lab[data.p$cond =="WT_DMSO"] <- "WT CTRL"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond ==cond2, paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

#cols = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
cols = c( "#E7298A","#c0c0c0", "#66A61E","#c0c0c0", "#E6AB02","#c0c0c0", "#7570B3", "#c0c0c0")

##  Plot domains
p <- ggplot(data.p, aes(x=PackEff, y=Rmax, col=lab2)) +
  scale_x_continuous(limits = c(0.3,1))+scale_y_continuous(limits = c(0,250))+
  geom_point( alpha=0.8, size = 3)+ theme_pubr()+
  scale_color_manual(values=cols, aesthetics = c("color", "fill"))+
  geom_vline(xintercept=pe, col="darkred", lty = "dashed") +theme( legend.position="left",plot.title = element_text(size=16), plot.subtitle = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=rad, col="darkred", lty = "dashed") + ylab("Domain Radius (nm)")+ xlab("Packing Efficiency")
p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.9, xfill = Rmax), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.9, yfill = PackEff), position = "dodge")+
  scale_xfill_manual(values =  cols) +
  scale_yfill_manual(values =  cols) + theme(ggside.panel.scale = .2)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 

p + geom_xsideboxplot(aes(y=lab2), varwidth=F, size     = 1 ,orientation = "y")  +
  geom_ysideboxplot(aes(x=lab2), varwidth=F, size     = 1 ,orientation = "x") +
  scale_xfill_manual(values =  cols) +
  scale_yfill_manual(values =  cols) + theme(ggside.panel.scale = .2)+
  scale_xsidey_discrete(guide = guide_axis(angle = 0)) + 
  scale_ysidex_discrete(guide = guide_axis(angle = 90)) +theme_ggside_void()+
  theme(legend.position="none", ggside.axis.text.x= element_blank(),
        ggside.axis.text.y= element_blank()
  ) 

## pol2.quadrant.pol2auxVSpol2dmso
## pol1.quadrant.pol1auxVSpol1dmso

##-------------------------------------------------------## Get barplot

library(ggside)

cond1 = "Pol1_6hrs_Aux"
cond2 = "WT_DMSO"

data.p <- data[data$cond ==cond1 | data$cond ==cond2 ,]
data.p$lab <- NA
pe <- median(data.p$PackEff[data.p$cond ==cond2])
rad <- median(data.p$Rmax[data.p$cond ==cond2])

##-- reorder data if needed for plotting layers

# Reorder rows by 'Name' column alphabetically from Z to A
data.p.sorted <- data.p[order(data.p$cond, decreasing = TRUE), ]
data.p <- data.p.sorted 

##-- end reorder

## Split data up into quarter
data.p$lab[data.p$PackEff <= pe & as.numeric(data.p$Rmax) >= rad] <- "Decaying"
data.p$lab[data.p$PackEff >= pe & as.numeric(data.p$Rmax) >= rad] <- "Large Mature"
data.p$lab[data.p$PackEff >= pe & as.numeric(data.p$Rmax) <= rad] <- "Small Mature"
data.p$lab[data.p$PackEff <= pe & as.numeric(data.p$Rmax) <= rad] <- "Nascent"
#data.p$lab[data.p$cond =="WT_DMSO"] <- "WT CTRL"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond == cond2, paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

labs<- unique(data.p$lab2)
a <- nrow(data.p[data.p$lab2 == labs[1],])
b <- nrow(data.p[data.p$lab2 == labs[2],])
c <- nrow(data.p[data.p$lab2 == labs[3],])
d <- nrow(data.p[data.p$lab2 == labs[4],])
e <- nrow(data.p[data.p$lab2 == labs[5],])
f <- nrow(data.p[data.p$lab2 == labs[6],])
g <- nrow(data.p[data.p$lab2 == labs[7],])
h <- nrow(data.p[data.p$lab2 == labs[8],])

## Pol 1
barplot <- data.frame(labs, domains = rbind(a,b,c,d,e,f,g,h),value = c(4, 2, 8, 6, 1, 7, 5, 3))

## Pol 2
barplot <- data.frame(labs, domains = rbind(a,b,c,d,e,f,g,h),value = c(2, 4, 8, 6, 1, 3, 5, 7))

## Pol 1 aux vs pol 1 dmso
barplot <- data.frame(labs, domains = rbind(a,b,c,d,e,f,g,h),value = c(1, 5, 3, 8, 2, 4, 6, 8))

cols = c( "#E7298A","#c0c0c0", "#66A61E","#c0c0c0", "#E6AB02","#c0c0c0", "#7570B3", "#c0c0c0")

barplot %>% arrange(as.numeric(value)) %>%
ggplot(aes(x=reorder(labs, value), y = domains, fill=as.factor(labs)), color = labs ) + 
  geom_bar( stat = "identity", color = "black", size = 1.5) +  theme_pubr()+
  scale_color_manual(values=rep("black",7))+
  scale_fill_manual(values=cols)+ scale_y_continuous(breaks = seq(0,60,5), labels= seq(0,60,5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  #theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x = element_blank()) + 
  ylab("Domains/µm^3")+ xlab("") + geom_text(aes(label = domains), size = 8, vjust = 1.5, colour = "black")


#####################################################################################
# ChromSTEM Pol2 plotting                              
#####################################################################################

library(ggforce)
library(PupillometryR)
library(here)
library(tidyverse)
here()
# Set directory variables here
dir.data <- here("Data/Wing Shun/Pol1 Chrom STEM") |> str_remove("/Pol2 Submission/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "alldomains_081725.csv"), sep= ",", header= TRUE))

brewer.pal(6, "Dark2")
cols <- c("#7570B3", "#E6AB02")

## Function to plot chromSTEM data
plot_violin <- function(data, title, xlab, ylab, y.pos,col.y, col.x, cols, comps) {
  
  require(ggpubr)
  require(ggplot2)
  
  my_comparisons <- comps
  
  p <- ggplot(data, aes_string(x={{col.x}}, y = {{col.y}},  color={{col.x}}, fill={{col.x}})) +
    geom_flat_violin(color = NA,trim = F, position = position_nudge(x = .15))+ 
    stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(width = 0.25,notch=TRUE, notchwidth = 0.8, outlier.shape = NA)+ scale_color_manual(values=rep("black", length(unique(data$cond))+1))+ scale_fill_manual(values = cols) +
    ggtitle(paste0(title)) + xlab("") + ylab(ylab)+ scale_x_discrete(labels = xlab)+
    theme_classic() + theme(legend.position="",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0),
                            panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) + 
    stat_compare_means(#label = "p.format",
      label = "p.signif",
      comparisons = my_comparisons, 
      #method = "t.test", 
      ref.group = ".all.", ## for multiple comparisons
      label.y = y.pos
    ) 
  
  print(p)
  
}

data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="WT_DMSO" ,]
data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="Pol2_DMSO" ,]

## single comparisons
plot_violin(data=data.p, title = "Chromatin Packing Scaling", ylab = "D", xlab=c("Aux", "DMSO"), y.pos = 2.95, col.y = "D", col.x="cond", cols= cols, comps = list(c("Pol2_6hrs_Aux", "Pol2_DMSO" )))
plot_violin(data=data.p, title = "Domain Radius", ylab = "Domain Radius (nm)", xlab=c("Aux", "DMSO"), y.pos = 225, col.y = "Rmax", col.x="cond", cols= cols, list(c("Pol2_6hrs_Aux", "Pol2_DMSO" )))
plot_violin(data=data.p, title = "Domain CVC", ylab = "Chromatin Volume Concentration", xlab=c("Aux", "DMSO"), y.pos = 0.5, col.y = "Density", col.x="cond", cols= cols, comps = list(c("Pol2_6hrs_Aux", "Pol2_DMSO"  ))) 

cols <- c("#7570B3", "#E6AB02","#66A61E")

## Multiple comparisons
plot_violin(data=data, title = "Chromatin Packing Scaling", ylab = "D", xlab=c("Pol 1 Aux","Pol 1 DMSO",  "WT DMSO"), y.pos = c(2.9, 2.95, 3), col.y = "D", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" ))) 
plot_violin(data=data, title = "Domain Radius", ylab = "Domain Radius (nm)", xlab=c("Pol 1 Aux","Pol 1 DMSO", "WT DMSO"), y.pos = c(210, 233, 245), col.y = "Rmax", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" ))) 
plot_violin(data=data, title = "Domain CVC", ylab = "Chromatin Volume Concentration", xlab=c("Pol 1 Aux","Pol 1 DMSO", "WT DMSO"), y.pos = c(0.45,0.5,0.55), col.y = "Density", col.x="cond", cols= cols, comps = list(c( "Pol1_DMSO", "WT_DMSO"  ),c( "Pol1_6hrs_Aux","Pol1_DMSO" ), c("Pol1_6hrs_Aux", "WT_DMSO" )))

##-------------------------------------------------------## Plot quadrants

library(ggside)
library(ggpubr)

data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="Pol2_DMSO" ,]
data.p$lab <- NA
pe <- median(data.p$PackEff[data.p$cond =="Pol2_DMSO"])
rad <- median(data.p$Rmax[data.p$cond =="Pol2_DMSO"])

## Split data up into quarter
data.p$lab[data.p$PackEff <= pe & data.p$Rmax >= rad] <- "Decaying"
data.p$lab[data.p$PackEff >= pe & data.p$Rmax >= rad] <- "Large Mature"
data.p$lab[data.p$PackEff >= pe & data.p$Rmax <= rad] <- "Small Mature"
data.p$lab[data.p$PackEff <= pe & data.p$Rmax <= rad] <- "Nascent"
#data.p$lab[data.p$cond =="WT_DMSO"] <- "WT CTRL"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond =="Pol2_DMSO", paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

#cols = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
cols = c( "#E7298A","#c0c0c0", "#66A61E","#c0c0c0", "#E6AB02","#c0c0c0", "#7570B3","#c0c0c0")

##  Plot domains
p <- ggplot(data.p, aes(x=PackEff, y=Rmax, col=lab2)) +
  scale_x_continuous(limits = c(0.3,1))+scale_y_continuous(limits = c(0,250))+
  geom_point( alpha=0.8, size = 3)+ theme_pubr()+
  scale_color_manual(values=cols, aesthetics = c("color", "fill"))+
  geom_vline(xintercept=pe, col="darkred", lty = "dashed") +theme( legend.position="left",plot.title = element_text(size=16), plot.subtitle = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=rad, col="darkred", lty = "dashed") + ylab("Domain Radius (nm)")+ xlab("Packing Efficiency")
p

p + geom_xsidedensity(aes(y=after_stat(density), alpha=0.9, xfill = Rmax), position = "dodge") +
  geom_ysidedensity(aes(x=after_stat(density), alpha=0.9, yfill = PackEff), position = "dodge")+
  scale_xfill_manual(values =  cols) +
  scale_yfill_manual(values =  cols) + theme(ggside.panel.scale = .2)+
  scale_xsidey_continuous() + scale_ysidex_continuous() 

p + geom_xsideboxplot(aes(y=lab2), varwidth=F, size     = 1 ,orientation = "y")  +
  geom_ysideboxplot(aes(x=lab2), varwidth=F, size     = 1 ,orientation = "x") +
  scale_xfill_manual(values =  cols) +
  scale_yfill_manual(values =  cols) + theme(ggside.panel.scale = .2)+
  scale_xsidey_discrete(guide = guide_axis(angle = 0)) + 
  scale_ysidex_discrete(guide = guide_axis(angle = 90)) +theme_ggside_void()+
  theme(legend.position="none", ggside.axis.text.x= element_blank(),
        ggside.axis.text.y= element_blank()
  ) 

##-------------------------------------------------------## Get barplot

library(ggside)
library(dplyr)

data.p <- data[data$cond =="Pol2_6hrs_Aux" | data$cond =="Pol2_DMSO" ,]
data.p$lab <- NA
pe <- median(data.p$PackEff[data.p$cond =="Pol2_DMSO"])
rad <- median(data.p$Rmax[data.p$cond =="Pol2_DMSO"])

## Split data up into quarter
data.p$lab[data.p$PackEff <= pe & as.numeric(data.p$Rmax) >= rad] <- "Decaying"
data.p$lab[data.p$PackEff >= pe & as.numeric(data.p$Rmax) >= rad] <- "Large Mature"
data.p$lab[data.p$PackEff >= pe & as.numeric(data.p$Rmax) <= rad] <- "Small Mature"
data.p$lab[data.p$PackEff <= pe & as.numeric(data.p$Rmax) <= rad] <- "Nascent"
#data.p$lab[data.p$cond =="WT_DMSO"] <- "WT CTRL"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond =="Pol2_DMSO", paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

labs<- unique(data.p$lab2)
a <- nrow(data.p[data.p$lab2 == labs[1],])
b <- nrow(data.p[data.p$lab2 == labs[2],])
c <- nrow(data.p[data.p$lab2 == labs[3],])
d <- nrow(data.p[data.p$lab2 == labs[4],])
e <- nrow(data.p[data.p$lab2 == labs[5],])
f <- nrow(data.p[data.p$lab2 == labs[6],])
g <- nrow(data.p[data.p$lab2 == labs[7],])
h <- nrow(data.p[data.p$lab2 == labs[8],])

values <- c("Large Mature KO" , "Large Mature WT","Small Mature KO", "Small Mature WT",
            "Nascent KO", "Nascent WT", "Decaying KO", "Decaying WT")
barplot <- data.frame(labs, domains = rbind(a,b,c,d,e,f,g, h),value = values)

## sort data frame
sort_vars <- c("labs", "value")
barplot_ordered <- barplot %>%
  arrange(across(all_of(sort_vars)))

barplot_ordered$factor <- c(7,8,1,2,5,6,3,4)


cols = c( "#E7298A","#c0c0c0", "#66A61E","#c0c0c0", "#E6AB02","#c0c0c0",  "#7570B3","#c0c0c0")

barplot_ordered %>% #arrange(as.numeric(value)) %>%
  ggplot(aes(x=reorder(labs, factor), y = domains, fill=as.factor(labs)), color = labs ) + 
  geom_bar( stat = "identity", color = "black", size = 1.5) +  theme_pubr()+
  scale_color_manual(values=rep("black",7))+
  scale_fill_manual(values=cols)+ scale_y_continuous(breaks = seq(0,60,5), labels= seq(0,60,5))+
  #theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x = element_blank()) + 
  ylab("Domains/µm^3")+ xlab("") + 
  geom_text(aes(label = domains), size = 8, vjust = 1.5, colour = "black")

##-------------------------------------------------------## Get barplot with proportions and stat testing

library(dplyr)
library(ggplot2)
library(binom) # Required for binom.confint

labs<- unique(data.p$lab2)
a <- nrow(data.p[data.p$lab2 == labs[1],])
b <- nrow(data.p[data.p$lab2 == labs[2],])
c <- nrow(data.p[data.p$lab2 == labs[3],])
d <- nrow(data.p[data.p$lab2 == labs[4],])
e <- nrow(data.p[data.p$lab2 == labs[5],])
f <- nrow(data.p[data.p$lab2 == labs[6],])
g <- nrow(data.p[data.p$lab2 == labs[7],])
h <- nrow(data.p[data.p$lab2 == labs[8],])


group <- c("Large Mature" , "Small Mature", "Nascent","Decaying")
bar <- data.frame(group,wins = as.numeric(c(a,b,c,d)), total = as.numeric(c(e,f,g, h)) )

add.later <- bar[c(2),]
bar<- bar[-c(2),]

# 2. Calculate Wilson Intervals using dplyr
 barplot <- bar%>%
  group_by(group) %>%
  mutate(proportion = wins/total) %>% 
  mutate(lower=prop.test(wins, total, correct = FALSE)$conf.int[2]) %>% 
  mutate( higher = prop.test(wins, total, correct = FALSE)$conf.int[1])

## add in extra data
barplot <- rbind(barplot, c(add.later, proportion = 1.2272, lower = NA, higher = NA) )
barplot$factor <- c(1,3,4,2)
cols = c( "#E7298A", "#66A61E", "#E6AB02",  "#7570B3")

barplot %>% #arrange(as.numeric(value)) %>%
  ggplot(aes(x=reorder(group, factor), y = proportion, fill=as.factor(group))) + 
    theme_pubr()+
  geom_bar( stat = "identity", color = "black", linewidth = 1.5) +
  scale_color_manual(values=rep("black",7))+
  geom_errorbar(aes(ymin = lower, ymax = higher), width = 0.2, size = 1,color = "black") +
  scale_fill_manual(values=cols)+ scale_y_continuous(breaks = seq(0,1.25,0.25), labels= seq(0,1.25,0.25))+
  theme(legend.position="none", text = element_text(size = 20)) + 
  ylab("Percentage")+ xlab("") + 
  geom_label(
    aes(label = round(proportion, 3) ),
    fill = "white",
    vjust = 1.5,
    color = "black", # Text color
    label.padding = unit(0.5, "lines"), # Adjust padding as needed
    label.size = 1 # Remove the border
  ) +scale_x_discrete(labels = function(group) str_wrap(group, width = 5))

#####################################################################################
# Plot quadrants for CVC versus size                      
#####################################################################################


library(ggside)
library(ggpubr)

cond1 = "Pol2_6hrs_Aux"
cond2 = "Pol2_DMSO"

data.p <- data[data$cond ==cond1 | data$cond ==cond2 ,]
data.p$lab <- NA

##-- reorder data if needed for plotting layers

# Reorder rows by 'Name' column alphabetically from Z to A
data.p.sorted <- data.p[order(data.p$cond, decreasing = TRUE), ]
data.p <- data.p.sorted 

##-- end reorder

## Split data up into quarter
data.p$lab[data.p$Rmax <= 80] <- "Small"
data.p$lab[data.p$Rmax >= 80 & data.p$Rmax <= 120] <- "Average"
data.p$lab[data.p$Rmax >= 120] <- "Large"


md.small <- median(data.p$Density[data.p$cond ==cond2 & data.p$lab == "Small"])
md.ave <- median(data.p$Density[data.p$cond ==cond2 & data.p$lab == "Average"])
md.large <- median(data.p$Density[data.p$cond ==cond2 & data.p$lab == "Large"])

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond ==cond2, paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

#cols = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
cols = c( "#E7298A","#c0c0c0", "#66A61E","#c0c0c0", "#7570B3","#c0c0c0")

##  Plot domains
p <- ggplot(data.p, aes(x=Density, y=Rmax, col=lab2)) +
  #scale_x_continuous(limits = c(0.3,1))+scale_y_continuous(limits = c(0,250))+
  #geom_point( alpha=0.8, size = 3)+ 
  geom_point(data = subset(data.p, cond == cond2), alpha=0.8, size = 3) + # Group B points first
  geom_point(data = subset(data.p, cond == cond1), alpha=0.8, size = 3) + 
  theme_pubr()+
  scale_color_manual(values=cols, aesthetics = c("color", "fill"))+
  geom_segment(aes(x = md.small, y = 5, xend = md.small, yend = 75), lty = "dashed", color = "#41045AFF", linewidth = 0.2)+
  geom_segment(aes(x = md.ave, y = 85, xend = md.ave, yend = 115), lty = "dashed", color = "#41045AFF", linewidth = 0.2)+
  geom_segment(aes(x = md.large, y = 125, xend = md.large, yend = max(Rmax)-5), lty = "dashed", color = "#41045AFF", linewidth = 0.2)+
  theme( legend.position="left",plot.title = element_text(size=16), plot.subtitle = element_text(size=12), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=c(80, 120), col="darkred", lty = "dashed") + ylab("Domain Radius (nm)")+ xlab("CVC")+
  annotate("text", x = 0.4, y = 0, label = "0-80 nm", color = "black", size = 5,hjust = "right")+
  annotate("text", x = 0.4, y = 90, label = "80-120 nm", color = "black", size = 5,hjust = "right")+
  annotate("text", x = 0.4, y = 130, label = "> 120 nm", color = "black", size = 5,hjust = "right")

p + geom_xsideboxplot(aes(y=lab2), varwidth=F, size     = 1 ,orientation = "y")  +
  geom_ysideboxplot(aes(x=lab2), varwidth=F, size     = 1 ,orientation = "x") +
  scale_xfill_manual(values =  cols) +
  scale_yfill_manual(values =  cols) + theme(ggside.panel.scale = .2)+
  scale_xsidey_discrete(guide = guide_axis(angle = 0)) + 
  scale_ysidex_discrete(guide = guide_axis(angle = 90)) +theme_ggside_void()+
  theme(legend.position="none", ggside.axis.text.x= element_blank(),
        ggside.axis.text.y= element_blank()
  ) 

## pol2.quadrant.pol2auxVSpol2dmso
## pol1.quadrant.pol1auxVSpol1dmso

##-------------------------------------------------------## Get barplot

#####################################################################################
# ChromSTEM Pol2 CVC plotting                              
#####################################################################################

library(ggforce)
library(PupillometryR)
library(here)
library(tidyverse)
library(ggside)
library(ggpubr)


# Set directory variables here
dir.data <- here("Data/Wing Shun/Pol1 Chrom STEM") |> str_remove("/Pol2 Submission/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "alldomains_081725.csv"), sep= ",", header= TRUE))

cond1 = "Pol2_6hrs_Aux"
cond2 = "Pol2_DMSO"

data.p <- data[data$cond ==cond1 | data$cond ==cond2 ,]
data.p$lab <- NA

## Split data up into quarter
data.p$lab[data.p$Rmax <= 80] <- "Small"
data.p$lab[data.p$Rmax >= 80 & data.p$Rmax <= 120] <- "Average"
data.p$lab[data.p$Rmax >= 120] <- "Large"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <- ifelse(data.p$cond ==cond2, paste0(data.p$lab, " WT"), paste0(data.p$lab, " KO"))

## Get CVC dist plots
library(ggdist)

p <- ggplot(data.p, aes(x = cond, y = Density)) +
  geom_point(
    aes(color = cond),
    position = position_jitter(width = 0.1, seed = 1), # Add jitter for visibility
    alpha = 0.5,
    size = 1.5
  ) +
  geom_boxplot(
    width = 0.55, 
    outlier.shape = NA, # Hide outliers since we show raw data
    alpha = 0.7
  ) +
  stat_summary(
    fun = mean, 
    geom = "line",   # Adds a point for the mean
    shape = 18,       # Diamond shape for the mean point
    size = 3,         # Size of the mean point
    color = "black"   # Color of the mean point
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Density Distributions By Domain Size",
    y = "CVC",
    x = ""
  ) +
  theme_minimal() +
  scale_x_discrete(labels =c("Aux","DMSO"))+
  theme_classic() + theme(legend.position="",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0),
                          panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) 


## add statistics in manually
library(rstatix)
## calculate stats manually
stat.test <-  data.p %>%
  group_by(lab) %>%
  rstatix::t_test(Density~cond ) %>% # Perform test for mpg vs am within each vs group
  rstatix::add_significance() %>%
  mutate(y.position = 0.4)

# Convert the 'category' variable to a factor and set the desired level order
## plot facet with manual stats
p +facet_wrap(~factor(lab, c("Small","Average", "Large"))) +
  stat_pvalue_manual(stat.test, inherit.aes = FALSE,
                     label = "p.signif" ) # Position the brackets


#####################################################################################
# Nico K9 plotting                        
#####################################################################################

# Set directory variables here
dir.data <- here("new data/Nico/k9_analysis/") |> str_remove("/code")
list.files(dir.data)

data <- data.frame(read.table(file.path(dir.data, "K9_concatenated_data_Full.csv"), sep= ",", header= TRUE))

colnames(data)
## -- prepare data
data.p <- data[,c("Auxin_density", "WT_density","Auxin_size_category" ,"WT_size_category")]

data.1 <- data[,c("Auxin_density", "Auxin_size_category")]
data.2 <- data[,c( "WT_density","WT_size_category")]

data.1 <- data[,c("Auxin_density", "WT_density")] 
data.2 <- data[,c("Auxin_size_category" ,"WT_size_category")]

colnames(data.1) <- c("Density", "Size")
colnames(data.2) <- c("Density", "Size")

data.1$cond <- "aux"
data.2$cond <- "wt"

data.p <- rbind(data.1, data.2)


## Split data up into quarter
data.p$lab <- NA
data.p$lab[data.p$Size == 0 | data.p$Size == 1] <- "Small"
data.p$lab[data.p$Size == 2] <- "Average"
data.p$lab[data.p$Size == 3 | data.p$Size == 4] <- "Large"

data.p <- na.omit(data.p)
data.p$lab2 <- NA
data.p$lab2 <-  paste0(data.p$lab," ", data.p$cond)

library(ggdist)

p <- ggplot(data.p, aes(x = factor(cond, levels = c("wt", "aux") ), y = Density)) +
  geom_point(
    aes(color = cond),
    position = position_jitter(width = 0.1, seed = 1), # Add jitter for visibility
    alpha = 0.3,
    size = 1
  ) +
  stat_halfeye(
    aes(fill = factor(cond)),
    adjust = 0.5,
    width = 0.6,
    justification = -0.27, # Nudge the density to one side
    .width = 0,
    point_colour = NA
  )+
  geom_boxplot(
    width = 0.23, 
    outlier.shape = NA, # Hide outliers since we show raw data
    alpha = 0.7
  ) +
  stat_summary(
    fun = mean, 
    geom = "line",   # Adds a point for the mean
    shape = 18,       # Diamond shape for the mean point
    size = 3,         # Size of the mean point
    color = "black"   # Color of the mean point
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Density Distributions By K9 Core Size",
    y = "Density (1/nm²)",
    x = ""
  ) +
  theme_minimal() +
  scale_x_discrete(labels =c("DMSO","Aux"))+
  theme_classic() + theme(legend.position="",plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),text = element_text(size=14, family="Arial"), axis.text.x = element_text(angle = 0),
                          panel.border = element_rect(linetype = "solid", colour = "black", fill = "NA", size = 0.5)) 

p
## add statistics in manually
library(rstatix)
## calculate stats manually
stat.test <-  data.p %>%
  group_by(lab) %>%
  rstatix::t_test(Density~cond ) %>% # Perform test for mpg vs am within each vs group
  rstatix::add_significance() %>%
  mutate(y.position = c(0.005, 0.004, 0.0075) )

mean_by_group <- data.p %>%
  group_by(lab2) %>%
  summarise(mean_value = mean(Density))

# Convert the 'category' variable to a factor and set the desired level order
## plot facet with manual stats
p +facet_wrap(~factor(lab, c("Small","Average", "Large")),scales = "free") +
  stat_pvalue_manual(stat.test, inherit.aes = FALSE,
                     label = "p.signif" )


