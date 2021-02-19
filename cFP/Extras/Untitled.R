#############################################################################
### Comparison of experimental (cFP) and simulated (MGM, PGM) trajectories
### DS3 and DS4 - FIG. 5
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
### DS1/6/7/9 - FIG. 9
###################################

# Combine experimental (cFP) and model (simulated) trajectories
## DS1/6/7/9 
## FIG S9
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
  ggsave("FIG_S9A_REVISION.pdf", width = 5, height = 16)

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
  ggsave("FIG_S9B_REVISION.pdf", width = 5, height = 16)

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

DS1_6_7_9_plot + ggsave("FIG_S9C_REVISION.pdf", width = 5, height = 16)

###################################
### Plot different subset
### DS8 - FIG. 10E
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
  ggsave("FIG_S10E_top_REVISION.pdf", width = 4, height = 4)

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
  ggsave("FIG_S10E_bottom_REVISION.pdf", width = 4, height = 4)

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
  ggsave("FIG_6F_REVISION.pdf", width = 6, height = 4)

