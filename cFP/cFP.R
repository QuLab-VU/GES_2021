setwd("~/git/GES_2020/Joint_functions/")
source("SumSE.R") # Summarize function to get the summary statistics
source("../DrugResponse/functionsDRC.R") # Several functions for drug-response data

setwd('../cFP/')
library(lubridate)
library(gplots)
library(diprate)
library(tidyxl)
library(readxl)
require(ggplot2)
require(Hmisc)
library(plyr)
library(ggpubr)
library(reshape2)

# ============================================================================================================
# ============================================================================================================

# Open all cFP data files
## All of these samples were seeded at 1 cell per well, grown for 1 week in RPMI 1640, 
## and 3um Erlotinib was added. The plates were imaged once a day for 1 week.
files = list.files(pattern = "M.csv")


###################################
### 48h normalization only
### See below for 72h normalization
###################################

CFPdata = data.frame()
CFPNL2 = data.frame()

# Make dataframe with all 384 well entries and samples from files
for(i in 1:length(files)){
  df = cellCountCV(read.csv(files[i], header=T))
  f = strsplit(files[i],"_")
  df$Sample = as.character(f[[1]][2])
  lwell = unique(df$Well)
  for(k in 1:length(lwell)){
    df[df$Well == lwell[k],"Seq"] = k 
  }
  # Get the log2 normalized counts for EACH well in a sample
  dfNL2 = compNL2(df, ntimepoint = 3)
  # dfNL2 = compNL2(df, ntimepoint = 4)
  # Remove outer wells - media volume was inconsistent
  dfNL2 = subset(dfNL2,!(dfNL2$Row%in%c("R01","R16")|dfNL2$Column%in%c("C01","C24")))
  # INF filter
  d <- dfNL2
  dnew <- do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))
  dfNL2 <- dnew[complete.cases(dnew),]
  # Filter out entries that had less than 50 cells at time of treatment
  ind <- which(dfNL2[dfNL2$Time==0,]$Count > 50)
  d_norm_i <- dfNL2[dfNL2$Time==0,]
  Well_keep <- d_norm_i[ind,]$Well
  dfNL2_n <- dfNL2[dfNL2$Well %in% Well_keep,]
  # Combine into a master dataframe
  # CFPdata = rbind(CFPdata,df)
  CFPNL2 = rbind(CFPNL2,dfNL2_n)
  rm(f,df,dfNL2,lwell)
}

# Pruning data for the times that were after the normalization -
# different for visualizations and DIP rate calculations 
## Added one more data point for DIP rate calculations

adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 40)
# adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 65)
adCFPNL2$Time = as.numeric(adCFPNL2$Time)
adCFPNL2$nl2 = as.numeric(adCFPNL2$nl2)
adCFPNL2$Seq = as.numeric(adCFPNL2$Seq)

# Remove outliers in DS7 - only for earlier normalization 
## DS7 had a slightly higher cell count before drug penetrance
# adCFPNL2_sub = subset(adCFPNL2, !(adCFPNL2$nl2 > 4 & adCFPNL2$Sample == "PC9-DS7"))

# Rename columns 
adCFPNL2_sub = adCFPNL2[,c("Count","Time","Well","Sample",
                           "Seq", "l2", "nl2")]

# Replace sample names with those in paper
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-parental", "PC9-VU")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS1", "DS1")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS3", "DS3")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS4", "DS4")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS6", "DS6")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS7", "DS7")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS8", "DS8")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS9", "DS9")

# Identify timespoints - used below and to subsample simulations
DS1_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time))[1]
DS3_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time))[1]
DS4_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time))[1]
DS6_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time))[1]
DS7_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time))[1]
DS8_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time))[1]
DS9_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time))[1]
VU_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-VU")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-VU")$Time))[1]
MGH_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-MGH")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-MGH")$Time))[1]
BR1_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-BR1")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-BR1")$Time))[1]

# Creating subsets of data
sublines <- c("DS1", "DS3", "DS4", "DS6", "DS7", "DS8", "DS9")
cFP_sublines <- subset(adCFPNL2_sub, Sample %in% sublines)
all_lines <- c("PC9-VU", "PC9-MGH", "PC9-BR1", "DS1", 
               "DS3", "DS4", "DS6", "DS7", "DS8", "DS9")
cFP_all <- subset(adCFPNL2_sub, Sample %in% all_lines)
## For use below
cFP_48hNorm <- cFP_all

# Creating separate dataframes for each sample in a common format
## variable = well id; nl2 = normalized log 2 cell count;
## Line = sample name (i.e., cell line, subline); 
## Type = Experimental (cFP) or Simulated (see more below)
DS1_exp <- subset(cFP_all, Sample == "DS1")
DS1_exp$Time <- rep(DS1_times, rep = length(unique(DS1_exp$Time)))
DS1_exp$variable = as.factor(rep(seq(length(unique(DS1_exp$Well))), each = length(unique(DS1_exp$Time))))
colnames(DS1_exp)[4] <- "Line"
DS1_exp$Type <- "Experimental"
DS1_exp <- DS1_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS3_exp <- subset(cFP_all, Sample == "DS3")
DS3_exp$Time <- rep(DS3_times, rep = length(unique(DS3_exp$Time)))
DS3_exp$variable = as.factor(rep(seq(length(unique(DS3_exp$Well))), each = length(unique(DS3_exp$Time))))
colnames(DS3_exp)[4] <- "Line"
DS3_exp$Type <- "Experimental"
DS3_exp <- DS3_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS4_exp <- subset(cFP_all, Sample == "DS4")
DS4_exp$Time <- rep(DS4_times, rep = length(unique(DS4_exp$Time)))
DS4_exp$variable = as.factor(rep(seq(length(unique(DS4_exp$Well))), each = length(unique(DS4_exp$Time))))
colnames(DS4_exp)[4] <- "Line"
DS4_exp$Type <- "Experimental"
DS4_exp <- DS4_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS6_exp <- subset(cFP_all, Sample == "DS6")
DS6_exp$Time <- rep(DS6_times, rep = length(unique(DS6_exp$Time)))
DS6_exp$variable = as.factor(rep(seq(length(unique(DS6_exp$Well))), each = length(unique(DS6_exp$Time))))
colnames(DS6_exp)[4] <- "Line"
DS6_exp$Type <- "Experimental"
DS6_exp <- DS6_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS7_exp <- subset(cFP_all, Sample == "DS7")
DS7_exp <- subset(DS7_exp, Seq != "326") # removed because not every time point captured
DS7_exp$Time <- rep(DS7_times, rep = length(unique(DS7_exp$Time)))
DS7_exp$variable = as.factor(rep(seq(length(unique(DS7_exp$Well))), each = length(unique(DS7_exp$Time))))
colnames(DS7_exp)[4] <- "Line"
DS7_exp$Type <- "Experimental"
DS7_exp <- DS7_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS8_exp <- subset(cFP_all, Sample == "DS8")
DS8_exp$Time <- rep(DS8_times, rep = length(unique(DS8_exp$Time)))
DS8_exp$variable = as.factor(rep(seq(length(unique(DS8_exp$Well))), each = length(unique(DS8_exp$Time))))
colnames(DS8_exp)[4] <- "Line"
DS8_exp$Type <- "Experimental"
DS8_exp <- DS8_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS9_exp <- subset(cFP_all, Sample == "DS9")
DS9_exp$Time <- rep(DS9_times, rep = length(unique(DS9_exp$Time)))
DS9_exp$variable = as.factor(rep(seq(length(unique(DS9_exp$Well))), each = length(unique(DS9_exp$Time))))
colnames(DS9_exp)[4] <- "Line"
DS9_exp$Type <- "Experimental"
DS9_exp <- DS9_exp[,c("Time", "variable", "nl2", "Line", "Type")]

VU_exp <- subset(cFP_all, Sample == "PC9-VU")
VU_exp$Time <- rep(VU_times, rep = length(unique(VU_exp$Time)))
VU_exp$variable = as.factor(rep(seq(length(unique(VU_exp$Well))), each = length(unique(VU_exp$Time))))
colnames(VU_exp)[4] <- "Line"
VU_exp$Type <- "Experimental"
VU_exp <- VU_exp[,c("Time", "variable", "nl2", "Line", "Type")]

MGH_exp <- subset(cFP_all, Sample == "PC9-MGH")
MGH_exp$Time <- rep(MGH_times, rep = length(unique(MGH_exp$Time)))
MGH_exp$variable = as.factor(rep(seq(length(unique(MGH_exp$Well))), each = length(unique(MGH_exp$Time))))
colnames(MGH_exp)[4] <- "Line"
MGH_exp$Type <- "Experimental"
MGH_exp <- MGH_exp[,c("Time", "variable", "nl2", "Line", "Type")]

BR1_exp <- subset(cFP_all, Sample == "PC9-BR1")
BR1_exp <- subset(BR1_exp, Seq != "271") # Removed because not every time point captured
BR1_exp$Time <- rep(BR1_times, rep = length(unique(BR1_exp$Time)))
BR1_exp$variable = as.factor(rep(seq(length(unique(BR1_exp$Well))), each = length(unique(BR1_exp$Time))))
colnames(BR1_exp)[4] <- "Line"
BR1_exp$Type <- "Experimental"
BR1_exp <- BR1_exp[,c("Time", "variable", "nl2", "Line", "Type")]

# Put cleaned data in common dataframe
cFP_all <- rbind(VU_exp, MGH_exp, BR1_exp, DS1_exp, DS3_exp,
                 DS4_exp, DS6_exp, DS7_exp, DS8_exp, DS9_exp)

# CFP proliferation trajectories for each PC9 cell line family member
## Colors and labels
cols_all <- c("DS1" = "coral", "DS3" = "brown", "DS4" = "deepskyblue", "DS6" = "deeppink",
              "DS7" = "darkorchid", "DS8" = "seagreen", "DS9" = "gold", "PC9-VU" = "blue",
              "PC9-MGH" = "green", "PC9-BR1" = "red")
pops_all <- c("DS1", "DS3", "DS4", "DS6", "DS7", 
              "DS8","DS9", "PC9-VU", "PC9-MGH", "PC9-BR1")
labs_all <- c("DS1", "DS3", "DS4", "DS6", "DS7", 
              "DS8","DS9", "PC9-VU", "PC9-MGH", "PC9-BR1")
levels(cFP_all$Line) <- c("PC9-VU", "PC9-MGH", "PC9-BR1", "DS1", 
                          "DS3", "DS4", "DS6", "DS7", "DS8", "DS9")

# Dataframe with the number of colonies (i.e., experimental replicated)
## Used to annotate plots
n_DF_all <- data.frame(Line = c("PC9-VU", "PC9-MGH", "PC9-BR1", "DS1", 
                                "DS3", "DS4", "DS6", "DS7", "DS8", "DS9"),
                       label=c(paste("n =", as.character(length(unique(subset(cFP_all, Line == "PC9-VU")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "PC9-MGH")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "PC9-BR1")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS1")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS3")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS4")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS6")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS7")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS8")$variable)))),
                               paste("n =", as.character(length(unique(subset(cFP_all, Line == "DS9")$variable))))))


###################################################
### END OF 48h normalization
### 72h normalization here
### Plots - FIG. 5A,5B,S3D,S3E,S4A,S9A,S9B,S10E
###################################################

CFPdata = data.frame()
CFPNL2 = data.frame()

# Make dataframe with all 384 well entries and samples from files
for(i in 1:length(files)){
  df = cellCountCV(read.csv(files[i], header=T))
  f = strsplit(files[i],"_")
  df$Sample = as.character(f[[1]][2])
  lwell = unique(df$Well)
  for(k in 1:length(lwell)){
    df[df$Well == lwell[k],"Seq"] = k 
  }
  # Get the log2 normalized counts for EACH well in a sample
  # dfNL2 = compNL2(df, ntimepoint = 3)
  dfNL2 = compNL2(df, ntimepoint = 4)
  # Remove outer wells - media volume was inconsistent
  dfNL2 = subset(dfNL2,!(dfNL2$Row%in%c("R01","R16")|dfNL2$Column%in%c("C01","C24")))
  # INF filter
  d <- dfNL2
  dnew <- do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))
  dfNL2 <- dnew[complete.cases(dnew),]
  # Filter out entries that had less than 50 cells at time of treatment
  ind <- which(dfNL2[dfNL2$Time==0,]$Count > 50)
  d_norm_i <- dfNL2[dfNL2$Time==0,]
  Well_keep <- d_norm_i[ind,]$Well
  dfNL2_n <- dfNL2[dfNL2$Well %in% Well_keep,]
  # Combine into a master dataframe
  # CFPdata = rbind(CFPdata,df)
  CFPNL2 = rbind(CFPNL2,dfNL2_n)
  rm(f,df,dfNL2,lwell)
}

# Pruning data for the times that were after the normalization -
# different for visualizations and DIP rate calculations 
## Added one more data point for DIP rate calculations

# adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 40)
adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 65)
adCFPNL2$Time = as.numeric(adCFPNL2$Time)
adCFPNL2$nl2 = as.numeric(adCFPNL2$nl2)
adCFPNL2$Seq = as.numeric(adCFPNL2$Seq)

# Remove outliers in DS7 - only for earlier normalization 
## DS7 had a slightly higher cell count before drug penetrance
# adCFPNL2_sub = subset(adCFPNL2, !(adCFPNL2$nl2 > 4 & adCFPNL2$Sample == "PC9-DS7"))

# Rename columns 
adCFPNL2_sub = adCFPNL2[,c("Count","Time","Well","Sample",
                           "Seq", "l2", "nl2")]

# Replace sample names with those in paper
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-parental", "PC9-VU")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS1", "DS1")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS3", "DS3")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS4", "DS4")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS6", "DS6")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS7", "DS7")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS8", "DS8")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS9", "DS9")

# Identify timespoints - used below and to subsample simulations
DS1_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time))[1]
DS3_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time))[1]
DS4_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time))[1]
DS6_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time))[1]
DS7_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time))[1]
DS8_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time))[1]
DS9_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time))[1]
VU_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-VU")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-VU")$Time))[1]
MGH_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-MGH")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-MGH")$Time))[1]
BR1_times <- round(unique(subset(adCFPNL2_sub, Sample == "PC9-BR1")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "PC9-BR1")$Time))[1]

# Creating subsets of data
sublines <- c("DS1", "DS3", "DS4", "DS6", "DS7", "DS8", "DS9")
cFP_sublines <- subset(adCFPNL2_sub, Sample %in% sublines)
all_lines <- c("PC9-VU", "PC9-MGH", "PC9-BR1", "DS1", 
               "DS3", "DS4", "DS6", "DS7", "DS8", "DS9")
cFP_all <- subset(adCFPNL2_sub, Sample %in% all_lines)

# Creating separate dataframes for each sample in a common format
## variable = well id; nl2 = normalized log 2 cell count;
## Line = sample name (i.e., cell line, subline); 
## Type = Experimental (cFP) or Simulated (see more below)
DS1_exp <- subset(cFP_all, Sample == "DS1")
DS1_exp$Time <- rep(DS1_times, rep = length(unique(DS1_exp$Time)))
DS1_exp$variable = as.factor(rep(seq(length(unique(DS1_exp$Well))), each = length(unique(DS1_exp$Time))))
colnames(DS1_exp)[4] <- "Line"
DS1_exp$Type <- "Experimental"
DS1_exp <- DS1_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS3_exp <- subset(cFP_all, Sample == "DS3")
DS3_exp$Time <- rep(DS3_times, rep = length(unique(DS3_exp$Time)))
DS3_exp$variable = as.factor(rep(seq(length(unique(DS3_exp$Well))), each = length(unique(DS3_exp$Time))))
colnames(DS3_exp)[4] <- "Line"
DS3_exp$Type <- "Experimental"
DS3_exp <- DS3_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS4_exp <- subset(cFP_all, Sample == "DS4")
DS4_exp$Time <- rep(DS4_times, rep = length(unique(DS4_exp$Time)))
DS4_exp$variable = as.factor(rep(seq(length(unique(DS4_exp$Well))), each = length(unique(DS4_exp$Time))))
colnames(DS4_exp)[4] <- "Line"
DS4_exp$Type <- "Experimental"
DS4_exp <- DS4_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS6_exp <- subset(cFP_all, Sample == "DS6")
DS6_exp$Time <- rep(DS6_times, rep = length(unique(DS6_exp$Time)))
DS6_exp$variable = as.factor(rep(seq(length(unique(DS6_exp$Well))), each = length(unique(DS6_exp$Time))))
colnames(DS6_exp)[4] <- "Line"
DS6_exp$Type <- "Experimental"
DS6_exp <- DS6_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS7_exp <- subset(cFP_all, Sample == "DS7")
DS7_exp <- subset(DS7_exp, Seq != "326") # removed because not every time point captured
DS7_exp$Time <- rep(DS7_times, rep = length(unique(DS7_exp$Time)))
DS7_exp$variable = as.factor(rep(seq(length(unique(DS7_exp$Well))), each = length(unique(DS7_exp$Time))))
colnames(DS7_exp)[4] <- "Line"
DS7_exp$Type <- "Experimental"
DS7_exp <- DS7_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS8_exp <- subset(cFP_all, Sample == "DS8")
DS8_exp$Time <- rep(DS8_times, rep = length(unique(DS8_exp$Time)))
DS8_exp$variable = as.factor(rep(seq(length(unique(DS8_exp$Well))), each = length(unique(DS8_exp$Time))))
colnames(DS8_exp)[4] <- "Line"
DS8_exp$Type <- "Experimental"
DS8_exp <- DS8_exp[,c("Time", "variable", "nl2", "Line", "Type")]

DS9_exp <- subset(cFP_all, Sample == "DS9")
DS9_exp$Time <- rep(DS9_times, rep = length(unique(DS9_exp$Time)))
DS9_exp$variable = as.factor(rep(seq(length(unique(DS9_exp$Well))), each = length(unique(DS9_exp$Time))))
colnames(DS9_exp)[4] <- "Line"
DS9_exp$Type <- "Experimental"
DS9_exp <- DS9_exp[,c("Time", "variable", "nl2", "Line", "Type")]

VU_exp <- subset(cFP_all, Sample == "PC9-VU")
VU_exp$Time <- rep(VU_times, rep = length(unique(VU_exp$Time)))
VU_exp$variable = as.factor(rep(seq(length(unique(VU_exp$Well))), each = length(unique(VU_exp$Time))))
colnames(VU_exp)[4] <- "Line"
VU_exp$Type <- "Experimental"
VU_exp <- VU_exp[,c("Time", "variable", "nl2", "Line", "Type")]

MGH_exp <- subset(cFP_all, Sample == "PC9-MGH")
MGH_exp$Time <- rep(MGH_times, rep = length(unique(MGH_exp$Time)))
MGH_exp$variable = as.factor(rep(seq(length(unique(MGH_exp$Well))), each = length(unique(MGH_exp$Time))))
colnames(MGH_exp)[4] <- "Line"
MGH_exp$Type <- "Experimental"
MGH_exp <- MGH_exp[,c("Time", "variable", "nl2", "Line", "Type")]

BR1_exp <- subset(cFP_all, Sample == "PC9-BR1")
BR1_exp <- subset(BR1_exp, Seq != "271") # Removed because not every time point captured
BR1_exp$Time <- rep(BR1_times, rep = length(unique(BR1_exp$Time)))
BR1_exp$variable = as.factor(rep(seq(length(unique(BR1_exp$Well))), each = length(unique(BR1_exp$Time))))
colnames(BR1_exp)[4] <- "Line"
BR1_exp$Type <- "Experimental"
BR1_exp <- BR1_exp[,c("Time", "variable", "nl2", "Line", "Type")]

# Put cleaned data in common dataframe
cFP_all <- rbind(VU_exp, MGH_exp, BR1_exp, DS1_exp, DS3_exp,
                 DS4_exp, DS6_exp, DS7_exp, DS8_exp, DS9_exp)


# Plots for cell line version experimental trajectories (not included elsewhere)
## PC9-BR1
ggplot(data = BR1_exp, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data = subset(n_DF_all, Line == "PC9-BR1"),
            aes(x = 30, y = -2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = cols_all,
                     labels = labs_all) +
  ylim(-2.5, 3) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_S2A.pdf", width = 5, height = 4)

## PC9-MGH
ggplot(data = MGH_exp, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data = subset(n_DF_all, Line == "PC9-MGH"),
            aes(x = 30, y = -2.9, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = cols_all,
                     labels = labs_all) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_S2B.pdf", width = 5, height = 4)


# Plot of VU cFP trajectories (grey) with DS lines population response inside (colors) 
## Load population trajectories from DrugResponse.R script
load("PopD_trajectories.RData")

## Subset to only sublines, and match timepoint structure
sublines_noDS8 <- subset(erlPC9Stat3_renorm, Sample != c("MGH", "VU"))
sublines_noDS8$Time <- sublines_noDS8$Time - unique(sublines_noDS8$Time)[1]
## This removes DS8 (reaches confluence quickly) and puts the sublines
## on the same timepoint structure as PC9-VU cFP data
sublines_noDS8 <- subset(sublines_noDS8, Time < 100)
## Create DS8 dataframe on same timepoint structure
subline_DS8 <- subset(erlPC9Stat3, Sample == "DS8")
subline_DS8 <- subset(subline_DS8, Time < 100)

## Changing color scheme so that PC9-VU is grey
cols_t <- c("DS1" = "coral", "DS3" = "brown", "DS4" = "deepskyblue", "DS6" = "deeppink",
            "DS7" = "darkorchid", "DS8" = "seagreen", "DS9" = "gold", "VU" = "grey50")
pops_t <- c("DS1", "DS3", "DS4", "DS6", "DS7", 
            "DS8","DS9", "VU")
labs_t <- c("DS1", "DS3", "DS4", "DS6", "DS7", 
            "DS8","DS9", "VU")

ggplot(VU_exp, aes(x=Time, y=nl2, group=variable)) +
  geom_line(size = 0.5, alpha = 0.25, colour = "grey50") +
  geom_line(data = sublines_noDS8, aes(x=Time, y=nl2, colour=Sample, group = Sample), size = 1) +
  geom_line(data = subline_DS8, aes(x=Time, y=nl2, colour=Sample, group = Sample), size = 1) +
  labs(x = "Time (hours) Post Drug Penetration", y = "Normalized Log2 Cell Count") +
  theme_bw() +
  geom_text(data = subset(n_DF_all, Line == "PC9-VU"),
            aes(x = 20, y = -2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_fill_manual(values = cols_t,
                    breaks = pops_t,
                    labels = labs_t) +
  scale_color_manual(values = cols_t,
                     breaks = pops_t,
                     labels = labs_t) +
  ylim(-2.5,3) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    legend.position = "right",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_S2C.pdf", width = 6.5, height = 4.5)


#############################################################################
### Comparison of experimental (cFP) and simulated (MGM, PGM) trajectories
### DS3 and DS4 - FIG. 5A-C
### Use 48 hour normalization
#############################################################################

# Pulling and cleaning model simulation data
setwd("~/git/GES_2020/cFP")
DS1_trajs <- read.csv('trajectories_DS1_G50_REVISION.csv', row.names = 1)
DS1_trajs <- DS1_trajs[,colSums(is.na(DS1_trajs))<nrow(DS1_trajs)]
names(DS1_trajs) <- seq(ncol(DS1_trajs))
DS1_trajs$Time <- seq(nrow(DS1_trajs)) -1
DS1_trajs <- melt(DS1_trajs, id = c("Time"), value.name = "nl2")
DS1_trajs$Line <- "DS1"
DS1_trajs$Type <- "Simulation"
DS1_trajs <- subset(DS1_trajs, Time %in% DS1_times)

DS3_trajs <- read.csv('trajectories_DS3_G50_REVISION.csv', row.names = 1)
DS3_trajs <- DS3_trajs[,colSums(is.na(DS3_trajs))<nrow(DS3_trajs)]
names(DS3_trajs) <- seq(ncol(DS3_trajs))
DS3_trajs$Time <- seq(nrow(DS3_trajs)) -1
DS3_trajs <- melt(DS3_trajs, id = c("Time"), value.name = "nl2")
DS3_trajs$Line <- "DS3"
DS3_trajs$Type <- "Simulation"
DS3_trajs <- subset(DS3_trajs, Time %in% DS3_times)

DS4_trajs <- read.csv('trajectories_DS4_G50_REVISION.csv', row.names = 1)
DS4_trajs <- DS4_trajs[,colSums(is.na(DS4_trajs))<nrow(DS4_trajs)]
names(DS4_trajs) <- seq(ncol(DS4_trajs))
DS4_trajs$Time <- seq(nrow(DS4_trajs)) -1
DS4_trajs <- melt(DS4_trajs, id = c("Time"), value.name = "nl2")
DS4_trajs$Line <- "DS4"
DS4_trajs$Type <- "Simulation"
DS4_trajs <- subset(DS4_trajs, Time %in% DS4_times)

DS6_trajs <- read.csv('trajectories_DS6_G50_REVISION.csv', row.names = 1)
DS6_trajs <- DS6_trajs[,colSums(is.na(DS6_trajs))<nrow(DS6_trajs)]
names(DS6_trajs) <- seq(ncol(DS6_trajs))
DS6_trajs$Time <- seq(nrow(DS6_trajs)) -1
DS6_trajs <- melt(DS6_trajs, id = c("Time"), value.name = "nl2")
DS6_trajs$Line <- "DS6"
DS6_trajs$Type <- "Simulation"
DS6_trajs <- subset(DS6_trajs, Time %in% DS6_times)

DS7_trajs <- read.csv('trajectories_DS7_G50_REVISION.csv', row.names = 1)
DS7_trajs <- DS7_trajs[,colSums(is.na(DS7_trajs))<nrow(DS7_trajs)]
names(DS7_trajs) <- seq(ncol(DS7_trajs))
DS7_trajs$Time <- seq(nrow(DS7_trajs)) -1
DS7_trajs <- melt(DS7_trajs, id = c("Time"), value.name = "nl2")
DS7_trajs$Line <- "DS7"
DS7_trajs$Type <- "Simulation"
DS7_trajs <- subset(DS7_trajs, Time %in% DS7_times)

DS8_trajs <- read.csv('trajectories_DS8_G50_REVISION.csv', row.names = 1)
DS8_trajs <- DS8_trajs[,colSums(is.na(DS8_trajs))<nrow(DS8_trajs)]
names(DS8_trajs) <- seq(ncol(DS8_trajs))
DS8_trajs$Time <- seq(nrow(DS8_trajs)) -1
DS8_trajs <- melt(DS8_trajs, id = c("Time"), value.name = "nl2")
DS8_trajs$Line <- "DS8"
DS8_trajs$Type <- "Simulation"
DS8_trajs <- subset(DS8_trajs, Time %in% DS8_times)

DS9_trajs <- read.csv('trajectories_DS9_G50_REVISION.csv', row.names = 1)
DS9_trajs <- DS9_trajs[,colSums(is.na(DS9_trajs))<nrow(DS9_trajs)]
names(DS9_trajs) <- seq(ncol(DS9_trajs))
DS9_trajs$Time <- seq(nrow(DS9_trajs)) -1
DS9_trajs <- melt(DS9_trajs, id = c("Time"), value.name = "nl2")
DS9_trajs$Line <- "DS9"
DS9_trajs$Type <- "Simulation"
DS9_trajs <- subset(DS9_trajs, Time %in% DS9_times)

# Calculate DIP rates for experimental data (cFP - 48h normalized)
library(lme4)
#For each sample, make a list of linear model fits at the normalized time point.
#Population Doublings per Hour of treatment (normalized to ~48 hours)
samples = unique(cFP_48hNorm$Sample)
DIPfits = data.frame()
for(i in 1:length(samples)){
  samp = subset(cFP_48hNorm,cFP_48hNorm$Sample == samples[i])
  fit = lmList(nl2 ~ Time | Well, data = samp)
  dfit = coef(fit)
  dfit$Well = row.names(dfit)
  row.names(dfit) = NULL
  dfit$Sample = samples[i]
  DIPfits = rbind(dfit, DIPfits)
}

DIPfits <- melt(DIPfits, id.vars = "Sample", measure.vars = "Time")
DIPfits <- DIPfits[,c("Sample", "value")]
names(DIPfits) <- c("Line", "DIP")

# Creating subsets of data
## DS3/4 - FIG. 5, DS1/6/7/9 - FIG. S9, DS8 - FIG. S10
DS3_4_data <- rbind(DS3_trajs, DS4_trajs, DS3_exp, DS4_exp)
DS1_6_7_9_data <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs,
                        DS1_exp, DS6_exp, DS7_exp, DS9_exp)

# Creating color-label subsets
cols_DS3_4 <- c("DS3" = "brown", "DS4" = "deepskyblue")
labs_DS3_4 <- c("DS3", "DS4")

# Combining experiments and simulations into common dataframe
DS3_4_sim <- rbind(DS3_trajs, DS4_trajs)
DS3_4_exp <- rbind(DS3_exp, DS4_exp)

# Plot experimental trajectories (DS3, DS4)
ggplot(data = DS3_4_exp, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  facet_wrap(.~Line, ncol = 1) +
  geom_text(data = subset(n_DF_all, Line %in% c("DS3", "DS4")),
            aes(x = 100, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey", "grey")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_5A_REVISION.pdf", width = 5, height = 8)

# Create model simulation division and death rate dataframe
## Annotations in plot
kdiv_DF_DS3_4 <- data.frame(Line = c("DS3", "DS4"),
                            label=c(paste("k[division] == 0.018"),
                                    paste("k[division] == 0.021"))) 
kdth_DF_DS3_4 <- data.frame(Line = c("DS3", "DS4"),
                            label=c(paste("k[death] == 0.0185"),
                                    paste("k[death] == 0.0178")))

# Plot simulated trajectories (DS3, DS4)
ggplot(data = DS3_4_sim, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv_DF_DS3_4,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth_DF_DS3_4,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("brown", "deepskyblue")) +
  facet_wrap(.~Line, ncol = 1) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_5B_REVISION.pdf", width = 5, height = 8)

# Load simulated distribution data (calculated from trajectories ploted above)
DS3_4_distribution <- read.csv("distributions_G50_REVISION.csv", row.names = 1)[,c(2,3)]
DS3_4_distribution <- melt(DS3_4_distribution)
names(DS3_4_distribution) <- c("Line", "DIP")
DS3_4_distribution$Type = "Simulated"

# Subset experimental DIP rates
DS3_4_distribution_exp <- subset(DIPfits, Line %in% c("DS3", "DS4"))
DS3_4_distribution_exp$Type = "Experimental"

# Combine dataframes
DS3_4_distribution_data <- rbind(DS3_4_distribution, DS3_4_distribution_exp)

# Load p-values (K-S test, see model code) for distribution comparison
## KS test
KS_pvalues <- read.csv("KSbootstrap_G50_REVISION.csv", row.names = 1)
pKS_DF_DS3_4 <- data.frame(Line = c("DS3", "DS4"),
                           label=c(paste("p =", as.character(round(mean(KS_pvalues$DS3),2)),
                                         "±", as.character(round(sd(KS_pvalues$DS3),2))),
                                   paste("p =", as.character(round(mean(KS_pvalues$DS4),2)),
                                         "±", as.character(round(sd(KS_pvalues$DS4),2)))))

## AD test
AD_pvalues <- read.csv("ADbootstrap_G50_REVISION.csv", row.names = 1)
pAD_DF_DS3_4 <- data.frame(Line = c("DS3", "DS4"),
                           label=c(paste("p =", as.character(round(mean(AD_pvalues$DS3),2)),
                                         "±", as.character(round(sd(AD_pvalues$DS3),2))),
                                   paste("p =", as.character(round(mean(AD_pvalues$DS4),2)),
                                         "±", as.character(round(sd(AD_pvalues$DS4),2)))))

DS3_4_plot <- ggplot() +
  facet_wrap(.~Line, ncol = 1) + xlim(-0.025, 0.020) +
  geom_density(data = DS3_4_distribution, aes(x=DIP, fill=Line, color=Line), 
               alpha=.25, size = 0.5) +  
  geom_density(data = DS3_4_distribution_exp, aes(x=DIP, fill=Line, color=Line), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("brown", "deepskyblue")) +
  scale_fill_manual(values = c("brown", "deepskyblue")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density") + 
  geom_text(data= pAD_DF_DS3_4,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) 

DS3_4_plot + ggsave("FIG_5C_REVISION.pdf", width = 5, height = 8)

###################################
### Plot different subset
### DS1/6/7/9 - FIG. S12
###################################

# Combine experimental (cFP) and model (simulated) trajectories
## DS1/6/7/9 
DS1_6_7_9_data <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs,
                        DS1_exp, DS6_exp, DS7_exp, DS9_exp)

# New color-line scheme
cols_DS1_6_7_9 <- c("DS1" = "coral", "DS6" = "deeppink",
                    "DS7" = "darkorchid", "DS9" = "gold")
labs_DS1_6_7_9 <- c("DS1", "DS6", "DS7", "DS9")

# Combine simulated and experimental subsets
DS1_6_7_9_sim <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs)
DS1_6_7_9_exp <- rbind(DS1_exp, DS6_exp, DS7_exp, DS9_exp)

# Plot subline experimental cFP data
ggplot(data = DS1_6_7_9_exp, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  facet_wrap(.~Line, ncol = 1) +
  geom_text(data = subset(n_DF_all, Line %in% c("DS1", "DS6", "DS7", "DS9")),
            aes(x = 100, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey", "grey", "grey", "grey")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_S12A_REVISION.pdf", width = 5, height = 16)

# Division and death rates associated with example simulations
## Annotated on model plots
kdiv_DF_DS1_6_7_9 <- data.frame(Line = c("DS1", "DS6", "DS7", "DS9"),
                                label=c(paste("k[division] == 0.017"),
                                        paste("k[division] == 0.028"),
                                        paste("k[division] == 0.016"),
                                        paste("k[division] == 0.014"))) 
kdth_DF_DS1_6_7_9 <- data.frame(Line = c("DS1", "DS6", "DS7", "DS9"),
                                label=c(paste("k[death] == 0.0156"),
                                        paste("k[death] == 0.0273"),
                                        paste("k[death] == 0.0141"),
                                        paste("k[death] == 0.0130")))

# Plot simulated data (colored by subline)
ggplot(data = DS1_6_7_9_sim, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv_DF_DS1_6_7_9,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth_DF_DS1_6_7_9,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  facet_wrap(.~Line, ncol = 1) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_S12B_REVISION.pdf", width = 5, height = 16)

# Load simulated distribution data (calculated from trajectories ploted above)
DS1_6_7_9_distribution <- read.csv("distributions_G50.csv", row.names = 1)[,c(1,4,5,6)]
DS1_6_7_9_distribution <- melt(DS1_6_7_9_distribution)
names(DS1_6_7_9_distribution) <- c("Line", "DIP")
DS1_6_7_9_distribution$Type = "Simulated"

# Subset experimental DIP rates
DS1_6_7_9_distribution_exp <- subset(DIPfits, Line %in% c("DS1", "DS6", "DS7", "DS9"))
DS1_6_7_9_distribution_exp$Type = "Experimental"

# Combine experimental (cFP) and simulated (model) distribution data
DS1_6_7_9_distribution_data <- rbind(DS1_6_7_9_distribution, DS1_6_7_9_distribution_exp)

## KS test
pKS_DF_DS1_6_7_9 <- data.frame(Line = c("DS1", "DS6", "DS7", "DS9"),
                               label=c(paste("p =", as.character(round(mean(KS_pvalues$DS1),2)),
                                             "±", as.character(round(sd(KS_pvalues$DS1),2))),
                                       paste("p =", as.character(round(mean(KS_pvalues$DS6),2)),
                                             "±", as.character(round(sd(KS_pvalues$DS6),2))),
                                       paste("p =", as.character(round(mean(KS_pvalues$DS7),2)),
                                             "±", as.character(round(sd(KS_pvalues$DS7),2))),
                                       paste("p =", as.character(round(mean(KS_pvalues$DS9),2)),
                                             "±", as.character(round(sd(KS_pvalues$DS9),2)))))

## AD test
pAD_DF_DS1_6_7_9 <- data.frame(Line = c("DS1", "DS6", "DS7", "DS9"),
                               label=c(paste("p =", as.character(round(mean(AD_pvalues$DS1),2)),
                                             "±", as.character(round(sd(AD_pvalues$DS1),2))),
                                       paste("p =", as.character(round(mean(AD_pvalues$DS6),2)),
                                             "±", as.character(round(sd(AD_pvalues$DS6),2))),
                                       paste("p =", as.character(round(mean(AD_pvalues$DS7),2)),
                                             "±", as.character(round(sd(AD_pvalues$DS7),2))),
                                       paste("p =", as.character(round(mean(AD_pvalues$DS9),2)),
                                             "±", as.character(round(sd(AD_pvalues$DS9),2)))))

# Plot experimental and simulated distribution comparisons
DS1_6_7_9_plot <- ggplot() +
  facet_wrap(.~Line, ncol = 1) +
  geom_density(data = DS1_6_7_9_distribution, aes(x=DIP, fill=Line, color=Line),
               alpha=.25, size = 0.5) +  
  geom_density(data = DS1_6_7_9_distribution_exp, aes(x=DIP), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  scale_fill_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density")+ xlim(-0.025, 0.020) +
  geom_text(data= pAD_DF_DS1_6_7_9,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14))

DS1_6_7_9_plot + ggsave("FIG_S12C_REVISION.pdf", width = 5, height = 16)

###################################
### Plot different subset
### DS8 - FIG. 5 (E-G)
###################################

# Plot experimental cFP trajectories for DS8
ggplot(data = DS8_exp, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data = subset(n_DF_all, Line == "DS8"),
            aes(x = 25, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey")) +
  ylim(-2.5,2.7) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_5E_REVISION.pdf", width = 4, height = 4)

# Add division and death rate information (PGM, 2 division and death rates) for DS8
kdiv1_DF <- data.frame(Line = c("DS8"),
                       label=c(paste("k[division] == 0.022"))) 
kdiv2_DF <- data.frame(Line = c("DS8"),
                       label=c(paste("k[division] == 0.024"))) 
kdth1_DF <- data.frame(Line = c("DS8"),
                       label=c(paste("k[death] == 0.0211"))) 
kdth2_DF <- data.frame(Line = c("DS8"),
                       label=c(paste("k[death] == 0.0169"))) 

# Plot DS8 (PGM) simulated trajectories
ggplot(data = DS8_trajs, aes(x=Time, y=nl2, group = variable, color = Line)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv1_DF,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth1_DF,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdiv2_DF,
            aes(x = 5,y = -1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth2_DF,
            aes(x = 5,y = -2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("seagreen")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_5F_REVISION.pdf", width = 4, height = 4)

# Load in simulated DS8 DIP rate distribution information
DS8_distribution <- read.csv("distributions_DS8_G50_REVISION.csv", row.names = 1)[,c(1)]
DS8_distribution <- as.data.frame(DS8_distribution)
DS8_distribution$Line <- "DS8"
DS8_distribution$Type = "Simulated"
names(DS8_distribution) <- c("DIP", "Line", "Type")

# Subset experimental DS8 experimental DIP rate distribution
DS8_distribution_exp <- subset(DIPfits, Line == "DS8")
DS8_distribution_exp$Type = "Experimental"

# Combine experimental and simulated data into common dataframe
DS8_distribution_data <- rbind(DS8_distribution, DS8_distribution_exp)

## KS test
KS_pvalues_DS8 <- read.csv("KSbootstrap_DS8_G50_REVISION.csv", row.names = 1)
pKS_DF_DS8 <- data.frame(Line = c("DS8"),
                         label=c(paste("p =", as.character(round(mean(KS_pvalues_DS8$DS8),2)),
                                       "±", as.character(round(sd(KS_pvalues_DS8$DS8),2)))))

## AD test
AD_pvalues_DS8 <- read.csv("ADbootstrap_DS8_G50_REVISION.csv", row.names = 1)
pAD_DF_DS8 <- data.frame(Line = c("DS8"),
                         label=c(paste("p =", as.character(round(mean(AD_pvalues_DS8$DS8),2)),
                                       "±", as.character(round(sd(AD_pvalues_DS8$DS8),2)))))


ggplot(DS8_distribution, aes(x=DIP, fill=Line, color=Line)) +
  geom_density(alpha=.25, size = 0.5) +  
  geom_density(data = DS8_distribution_exp, aes(x=DIP), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("seagreen")) +
  scale_fill_manual(values = c("seagreen")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density")+ xlim(-0.025, 0.020) + ylim(0,210) +
  # facet_grid(.~Sample) +
  geom_text(data= pAD_DF_DS8,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("FIG_5G_REVISION.pdf", width = 6, height = 4)


################################################################
### Plotting cohort DIP rate distribution comparisons
### Cell Line Versions, Sublines - FIGS. 2B and 2D
################################################################

# Establish cell line version (CLV) cohort
CLV = c("PC9-VU", "PC9-MGH", "PC9-BR1")

# Subset DIP rates to only CLV
CLV_DIP = subset(DIPfits, DIPfits$Line%in%CLV)

# CLV DIP rate distribution plotting
plt_CLV_cFP <- ggplot(CLV_DIP, aes(x=DIP, fill=Line, colour=Line)) +
  geom_density(alpha=.5, size = 1)+
  geom_vline(xintercept = 0, size = 1, colour = "black",
             linetype = "dashed") +
  theme_bw() + xlim(-0.03, 0.05) +
  labs(x = "DIP Rate", y = "Density")+
  scale_color_manual(values = c("red", "green", "blue"),
                     labels = c("PC9-BR1", "PC9-MGH", "PC9-VU")) +
  scale_fill_manual(values = c("red", "green", "blue"),
                    labels = c("PC9-BR1", "PC9-MGH", "PC9-VU")) +
  theme(
    legend.position = "none", legend.text = element_text(size=14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plt_CLV_cFP + ggsave("FIG_2B_revised.svg", width = 6, height = 5)
# # Change legend.position to right, then run arguments to get legend
# leg_CLV_cFP <- get_legend(plt_CLV_cFP)
# as_ggplot(leg_CLV_cFP) + ggsave("legend_FIG_2B.svg", width = 3, height = 3)

# Establish discrete subline (DS) cohort
sublines = c("DS1", "DS3", "DS4", "DS6",
             "DS7", "DS8", "DS9")

# Subset DIP rates to only sublines
subline_DIP = subset(DIPfits, DIPfits$Line%in%sublines)

# Subline DIP rate distribution plotting
plt_sublines_cFP <- ggplot(subline_DIP, aes(x=DIP, fill=Line, color=Line)) +
  geom_density(alpha=.1, size = 1) +
  geom_density(data = subset(CLV_DIP, Line == "PC9-VU"), aes(x=DIP), alpha = 0,
               size = 1, linetype = "dashed") +
  # geom_vline(xintercept = 0, size = 1, colour = "black",
  #            linetype = "dashed") +
  theme_bw() + xlim(-0.015, 0.025) +
  # facet_wrap(.~Line, ncol = 4) +
  labs(x = "DIP Rate", y = "Density")+
  scale_color_manual(values = c("coral", "brown", "deepskyblue",
                                "deeppink", "darkorchid", "seagreen",
                                "gold", "black")) +
  scale_fill_manual(values = c("coral", "brown", "deepskyblue",
                               "deeppink", "darkorchid", "seagreen",
                               "gold", "black")) +
  theme(
    legend.position = "none", legend.text = element_text(size=14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plt_sublines_cFP + ggsave("FIG_2D_revised.svg", width = 6, height = 5)
# # Change legend.position to right, then run arguments to get legend
# leg_sublines_cFP <- get_legend(plt_sublines_cFP)
# as_ggplot(leg_sublines_cFP) + ggsave("legend_FIG_2D.svg", width = 3, height = 3)


# Run statistical test (Mood's median test) to identify if subline distributions are distinct
library(RVAideMemoire)
pval_medians <- mood.medtest(DIP ~ Line,
                             data  = subline_DIP,
                             exact = FALSE)[5]