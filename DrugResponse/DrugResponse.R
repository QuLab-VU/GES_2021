setwd("~/git/GES_2020/Joint_functions/")
source("SumSE.R") # Summarize function to get the summary statistics;
source("utMA.R") # Time-averaging code

setwd('../DrugResponse/')
source("functionsDRC.R") # Several functions for drug-response data
library(diprate) # Library for DIP rate calculations (Nature Methods 13, 497–500(2016))

# Libraries for plotting and manipulating data
library(lubridate)
library(tidyxl)
library(readxl)
library(dplyr)
require(ggplot2)
require(Hmisc)
library(ggpubr)

# ============================================================================================================
# ============================================================================================================

# Read in untreated and Erlotinib treated datasets
untreatCVdata = read.csv("Parental-VU-MGH-BR1_Clones-DS1-3-4-6-7-8-9_DMSO.csv", header=T)
erlCVdata = read.csv("Parental-VU-MGH-BR1_Clones-DS1-3-4-6-7-8-9_Erl-3uM.csv", header=T)

# Convert pre-processed data to counts
untreatPC9 = cellCountCV(untreatCVdata,cellnumcol = "Cell.Nucleus", normPos = 1)
erlPC9 = cellCountCV(erlCVdata,cellnumcol = "Cell.Nucleus", normPos = 1)

# Upload platemap files - differentiate between sample and condition
untreatPfile="Platemaps/Untreat PC9 plate map.xlsx"
erlPfile="Platemaps/Erlotnib PC9 plate map.xlsx"

# Add labels to row-column entries
## Untreated
untreatGroups = getPlateGroupString96(untreatPfile)
untreatPC9L = styleDrugPlate96(untreatPfile, untreatPC9)

## Treated
erlGroups = getPlateGroupString96(erlPfile)
erlPC9L = styleDrugPlate96(erlPfile, erlPC9)

# Normalize for change in cell counts after drug change - Erlotinib treated
dcErlPC9L = remDrugChange(erlPC9L)

# Exponential smoothing of count data - Erlotinib treated
erlESmooth = expSmooth(dcErlPC9L, alpha = 0.6, countcol = "Count")

# Implementation of time-step moving average 
d = dcErlPC9L
a <- sapply(unique(d$Well), function(w) 
  round(utMA(d[d$Well==w,'Count'],d[d$Well==w,'Time'],tau=16),0))

# Add averaged count to dataframe
af = data.frame(a)
df = data.frame()
for(i in 1:length(unique(d$Well))){
  shub = subset(d,d$Well%in%unique(d$Well)[i])
  for(j in 1:length(unique(d$Time))){
    shub[shub$Time==unique(d$Time)[j], "sCount"] = af[j,i]
  }
  df = rbind(df,shub)
}


erlPC9SMA = df

# Compile statistics on for each timepoint-sample
## Untreated
unPC9Stat = compNL2Stats(untreatPC9L, untreatGroups, time=0)

## Untreated (non-smoothed)
erlPC9Stat = compNL2Stats(dcErlPC9L, erlGroups, time=0)

## Erlotinib treated (moving-average)
erlPC9Stat2 = compNL2Stats(erlPC9SMA, erlGroups, CountVar = "sCount", time=0)

## Erlotinib treated (exponential smoothing)
### Remove DS8/BR1 (to be added back - they reach confluence quickly b/c resistant to Erlotinib)
erlESmooth_cutout <- subset(erlESmooth, !erlESmooth$Sample %in% c("BR1","DS8"))
erlESmooth_BR1 <- subset(erlESmooth, Sample %in% "BR1")
erlESmooth_BR1 <- subset(erlESmooth_BR1, Time < 125)
erlESmooth_DS8 <- subset(erlESmooth, Sample %in% "DS8")
erlESmooth_DS8 <- subset(erlESmooth_DS8, Time < 125)
### Add DS8/BR1 back in (up to confluence time)
erlESmooth_cut <- rbind(erlESmooth_cutout, erlESmooth_BR1, erlESmooth_DS8)
erlPC9Stat3 = compNL2Stats(erlESmooth_cut, erlGroups, CountVar = "esCount", time = 0)

# Normalize data to time point after drug takes effect (~125h)
erlESmooth_renorm <- erlESmooth_cutout
erlPC9Stat3_renorm = compNL2Stats(erlESmooth_renorm, erlGroups, CountVar = "esCount", time = 125.1458)
erlPC9Stat3_renorm = subset(erlPC9Stat3_renorm, Time > 124)

save(erlPC9Stat3, erlPC9Stat3_renorm, file = "../cFP/PopD_trajectories.RData")

## Only plot from < 75h (before any wells reach confluence)
unPC9Stat_cut <- unPC9Stat[unPC9Stat$Time < 75,]

# Subset samples (for later plots)
smoothsub = subset(erlPC9Stat3, !erlPC9Stat3$Sample%in%c("BR1","MGH","VU"))
smoothsub_VU = subset(erlPC9Stat3, erlPC9Stat3$Sample%in%c("VU"))
smoothsub_CLV = subset(erlPC9Stat3, erlPC9Stat3$Sample%in%c("BR1","MGH","VU"))
untreated_CLV <- subset(unPC9Stat_cut, unPC9Stat_cut$Sample%in%c("BR1","MGH","VU"))
untreated_sublines <- subset(unPC9Stat_cut, !unPC9Stat_cut$Sample%in%c("BR1","MGH","VU"))

# Plot exponentially smoothed drug response for cell line versions
plt_CLV <- ggplot() +
  geom_point(data = smoothsub_CLV, aes(x=Time, y=nl2, color = Sample, group = Sample, fill = Sample)) +
  geom_smooth(data = smoothsub_CLV, aes(x=Time, y=nl2, color = Sample, group = Sample, fill = Sample)) +
  geom_smooth(data = untreated_CLV, aes(x=Time, y=nl2, color = Sample), 
              method = "lm", linetype = "dashed", se = FALSE) +
  theme_bw() + xlab("Time (hours)") + ylab("Normalized Log2 Cell Count") +
  scale_color_manual(values = c("red", "green", "blue")) +
  scale_fill_manual(values = c("red", "green", "blue")) +
  ylim(-1,4) +
  theme(legend.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text=element_text(size=14), legend.title = element_text(size=14),
        legend.position = "none", axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plt_CLV + ggsave("FIG_2A_revised.svg", width = 5, height = 5)
# # Change legend.position to right, then run arguments to get legend
# leg_CLV <- get_legend(plt_CLV)
# as_ggplot(leg_CLV) + ggsave("legend_FIG_2A.svg", width = 3, height = 1)

# Plot exponentially smoothed drug response for sublines
plt_sublines <- ggplot() +
  geom_point(data = smoothsub, aes(x=Time, y=nl2, color = Sample, group = Sample, fill = Sample)) +
  geom_smooth(data = smoothsub, aes(x=Time, y=nl2, color = Sample, group = Sample, fill = Sample)) +
  geom_smooth(data = subset(smoothsub_CLV, Sample == "VU"), 
              aes(x=Time, y=nl2, color = Sample), se = FALSE) +
  geom_smooth(data = untreated_sublines, aes(x=Time, y=nl2, color = Sample), 
              method = "lm", linetype = "dashed", se = FALSE) +
  theme_bw() + 
  xlab("Time (hours)") + ylab("Normalized Log2 Cell Count") + 
  # facet_wrap(Drug~Sample, ncol = 4) + 
  ylim(-1,4) +
  scale_color_manual(values = c("coral", "brown", "deepskyblue",
                                "deeppink", "darkorchid", "seagreen",
                                "gold", "black")) +
  scale_fill_manual(values = c("coral", "brown", "deepskyblue",
                               "deeppink", "darkorchid", "seagreen",
                               "gold", "black")) +
  theme(legend.text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14),
        legend.title = element_text(size=14),
        legend.position = "none",
        axis.title=element_text(size=14),
        # legend.position = c(0.95, -0.05), legend.justification = c(1, 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plt_sublines + ggsave("FIG_2C_revised.svg", width = 5, height = 5)
# Change legend.position to right, then run arguments to get legend
leg_sublines <- get_legend(plt_sublines)
as_ggplot(leg_sublines) + ggsave("legend_FIG_2C.svg", width = 3, height = 3)