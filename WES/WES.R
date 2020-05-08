setwd('~/Documents/QuarantaLab/GES_2020/WES/')

# devtools::install_github(repo="knausb/vcfR")
library(reshape2)
library(ggplot2)
library(devtools)
library(vcfR)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(UpSetR)
library(scales)
library(waffle)
library(purrr)
library(stringr)
library(ggrepel)
library(forcats)
library(scales)
library(GenVisR)
library(ggbio)

# Numbers obtained from VCFtools and compiled into common dataset
num_muts <- read.csv("number_mutations.csv", header = T)
num_muts$Sample <- c("VU", "MGH", "BR1", "DS3", 
                     "DS6", "DS7", "DS8", "DS9")

# Put data in format for ggplot 
num_muts <- melt(num_muts, id.vars = "Sample")
num_muts = num_muts[!(num_muts$variable == "Total"),]

# Plot total mutation (breakdown by SNPs/InDels) counts
ggplot(num_muts, aes(x=Sample, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("grey40", "grey90"), labels = c("SNPs", "InDels")) +
  theme_bw() + xlab("Population") + ylab("Number of Mutations") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_S5A.pdf", width = 12, height = 4)

# Initiate objects for data quality control (QC) analysis 
vcf <- read.vcfR("/Users/Corey/Documents/QuarantaLab/PC9/WXS/accre_vars/samples_called_vars_named.vcf.gz", verbose = TRUE)
dna <- ape::read.dna("/Volumes/quaranta/Data/WXS/bwa_ref_genome/Homo_sapiens_assembly38.fasta", format = "fasta")
gff <- read.table("/Volumes/quaranta/Data/RNAseq/PC9_scRNAseq/genes.gtf", sep="\t", quote="")

# Plot QC metrics for compiled dataset (all variants in all samples) - FIG. S5B
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)
chromoqc(chrom)

# Plot unique and intersections of mutations for each cohort
## from vcftools vcf-compare data: compiled in csvs provided
### Cell Line Versions (CLV)
shared_vars_CLV <- read.csv("shared_variants_CLV.csv", header = T)
upset_CLV <- c("VU" = shared_vars_CLV[1, "Number"],
               "MGH" = shared_vars_CLV[2, "Number"],
               "BR1" = shared_vars_CLV[3, "Number"],
               "VU&MGH" = shared_vars_CLV[4, "Number"],
               "VU&BR1" = shared_vars_CLV[5, "Number"],
               "MGH&BR1" = shared_vars_CLV[6, "Number"],
               "VU&MGH&BR1" = shared_vars_CLV[7, "Number"])

### FIG. S5C - saved as 6x9 PDF
upset(fromExpression(upset_CLV), order.by = "freq",
      sets.bar.color = c("green", "red", "blue"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 101000,
      mainbar.y.label = "Number of Mutations", sets.x.label = "Set Size", set_size.scale_max = 160000,
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

### Sublines excluding DS8
shared_vars_sublines_noDS8 <- read.csv("shared_variants_sublines_noDS8.csv", header = T)
upset_sublines_noDS8 <- c("DS3" = shared_vars_sublines_noDS8[1, "Number"],
                          "DS6" = shared_vars_sublines_noDS8[2, "Number"],
                          "DS7" = shared_vars_sublines_noDS8[3, "Number"],
                          "DS9" = shared_vars_sublines_noDS8[4, "Number"],
                          "DS3&DS6" = shared_vars_sublines_noDS8[5, "Number"],
                          "DS3&DS7" = shared_vars_sublines_noDS8[6, "Number"],
                          "DS3&DS9" = shared_vars_sublines_noDS8[7, "Number"],
                          "DS6&DS7" = shared_vars_sublines_noDS8[8, "Number"],
                          "DS6&DS9" = shared_vars_sublines_noDS8[9, "Number"],
                          "DS7&DS9" = shared_vars_sublines_noDS8[10, "Number"],
                          "DS3&DS6&DS7" = shared_vars_sublines_noDS8[11, "Number"],
                          "DS3&DS6&DS9" = shared_vars_sublines_noDS8[12, "Number"],
                          "DS3&DS7&DS9" = shared_vars_sublines_noDS8[13, "Number"],
                          "DS6&DS7&DS9" = shared_vars_sublines_noDS8[14, "Number"],
                          "DS3&DS6&DS7&DS9" = shared_vars_sublines_noDS8[15, "Number"])

### FIG. S5D - saved as 6x9 PDF
upset(fromExpression(upset_sublines_noDS8), order.by = "freq",
      sets.bar.color = c("brown", "deeppink", "gold", "darkorchid"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 106000,
      mainbar.y.label = "Number of Mutations", sets.x.label = "Set Size", set_size.scale_max = 160000,
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

### Sublines including DS8
shared_vars_sublines <- read.csv("shared_variants_sublines.csv", header = T)
upset_sublines <- c("DS3" = shared_vars_sublines[1, "Number"],
                    "DS6" = shared_vars_sublines[2, "Number"],
                    "DS7" = shared_vars_sublines[3, "Number"],
                    "DS8" = shared_vars_sublines[4, "Number"],
                    "DS9" = shared_vars_sublines[5, "Number"],
                    "DS3&DS6" = shared_vars_sublines[6, "Number"],
                    "DS3&DS7" = shared_vars_sublines[7, "Number"],
                    "DS3&DS8" = shared_vars_sublines[8, "Number"],
                    "DS3&DS9" = shared_vars_sublines[9, "Number"],
                    "DS6&DS7" = shared_vars_sublines[10, "Number"],
                    "DS6&DS8" = shared_vars_sublines[11, "Number"],
                    "DS6&DS9" = shared_vars_sublines[12, "Number"],
                    "DS7&DS8" = shared_vars_sublines[13, "Number"],
                    "DS7&DS9" = shared_vars_sublines[14, "Number"],
                    "DS8&DS9" = shared_vars_sublines[15, "Number"],
                    "DS3&DS6&DS7" = shared_vars_sublines[16, "Number"],
                    "DS3&DS6&DS8" = shared_vars_sublines[17, "Number"],
                    "DS3&DS6&DS9" = shared_vars_sublines[18, "Number"],
                    "DS3&DS7&DS8" = shared_vars_sublines[19, "Number"],
                    "DS3&DS7&DS9" = shared_vars_sublines[20, "Number"],
                    "DS3&DS8&DS9" = shared_vars_sublines[21, "Number"],
                    "DS6&DS7&DS8" = shared_vars_sublines[22, "Number"],
                    "DS6&DS7&DS9" = shared_vars_sublines[23, "Number"],
                    "DS6&DS8&DS9" = shared_vars_sublines[24, "Number"],
                    "DS7&DS8&DS9" = shared_vars_sublines[25, "Number"],
                    "DS3&DS6&DS7&DS8" = shared_vars_sublines[26, "Number"],
                    "DS3&DS6&DS7&DS9" = shared_vars_sublines[27, "Number"],
                    "DS3&DS6&DS8&DS9" = shared_vars_sublines[28, "Number"],
                    "DS3&DS7&DS8&DS9" = shared_vars_sublines[29, "Number"],
                    "DS6&DS7&DS8&DS9" = shared_vars_sublines[30, "Number"],
                    "DS3&DS6&DS7&DS8&DS9" = shared_vars_sublines[31, "Number"])

### FIG. S5E - saved as 6x9 PDF
upset(fromExpression(upset_sublines), order.by = "freq",
      sets.bar.color = c("brown", "seagreen", "deeppink", "gold", "darkorchid"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 106000,
      mainbar.y.label = "Number of Mutations", sets.x.label = "Set Size", set_size.scale_max = 160000,
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

# Circle plot of mean-centered mutation counts
## Clear other packages before running this code
library(tidyverse)
library(reshape2)
## Open and rename mutation count by chromosome
### Pulled from VEP summary HTML files and compiled into csv
muts_Chrom <- read.csv("mutations_byChromosome.csv", header = T)
colnames(muts_Chrom) <- c("Chromosome", "VU", "MGH", "BR1",
                          "DS3", "DS6", "DS7", "DS8", "DS9")
muts_Chrom_named <- muts_Chrom[,-1]
rownames(muts_Chrom_named) <- muts_Chrom[,1]

### Function to mean-center mutation counts
center_apply <- function(x) {
  apply(x, 1, function(y) y - mean(y))
}

### Mean center mutation count across cohort (CLV, sublines) by chromosome
#### CLV
muts_Chrom_named_center_CLV <- t(center_apply(muts_Chrom_named[,c("VU", "MGH", "BR1")]))
muts_Chrom_named_center_melt_CLV <- melt(data = muts_Chrom_named_center_CLV,
                                         id.vars = rownames(muts_Chrom_named_center_CLV),
                                         measure.vars = c("VU", "MGH", "BR1"))
#### Sublines excluding DS8
muts_Chrom_named_center_sublines_noDS8 <- t(center_apply(muts_Chrom_named[,c("DS3", "DS6", "DS7", "DS9")]))
muts_Chrom_named_center_melt_sublines_noDS8 <- melt(data = muts_Chrom_named_center_sublines_noDS8,
                                              id.vars = rownames(muts_Chrom_named_center_sublines_noDS8),
                                              measure.vars = c("DS3", "DS6", "DS7", "DS8", "DS9"))
#### Sublines including DS8
muts_Chrom_named_center_sublines <- t(center_apply(muts_Chrom_named[,c("DS3", "DS6", "DS7", "DS8", "DS9")]))
muts_Chrom_named_center_melt_sublines <- melt(data = muts_Chrom_named_center_sublines,
                                            id.vars = rownames(muts_Chrom_named_center_sublines),
                                            measure.vars = c("DS3", "DS6", "DS7", "DS8", "DS9"))

### Creating the circle plot - CLV
empty_bar <- 3 # Used for spacing
#### Set a number of 'empty bar' to add at the end of each group
data_CLV <- muts_Chrom_named_center_melt_CLV
colnames(data_CLV) <- c("group", "population", "value")
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data_CLV$group), ncol(data_CLV)) )

#### Add identifier label to dataframe
colnames(to_add) <- colnames(data_CLV)
to_add$group <- rep(levels(data_CLV$group), each=empty_bar)
data_CLV <- rbind(data_CLV, to_add)
data_CLV <- data_CLV %>% arrange(group)
data_CLV$id <- seq(1, nrow(data_CLV))

# Get the name and the y position of each label
label_data_CLV <- data_CLV
number_of_bar <- nrow(label_data_CLV)
angle <- 90 - 360 * (label_data_CLV$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data_CLV$hjust <- ifelse( angle < -90, 1, 0)
label_data_CLV$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a dataframe for base lines
base_data_CLV <- data_CLV %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a dataframe for grid (scales)
grid_data_CLV <- base_data_CLV
grid_data_CLV$end <- grid_data_CLV$end[ c( nrow(grid_data_CLV), 1:nrow(grid_data_CLV)-1)] + 1
grid_data_CLV$start <- grid_data_CLV$start - 1
grid_data_CLV <- grid_data_CLV[-1,]

ggplot(data_CLV, aes(x=as.factor(id), y=value, fill=population)) + 
  geom_bar(aes(x=as.factor(id), y=value, fill=population), stat="identity", alpha=0.5) +
  # Add reference lines
  geom_segment(data=grid_data_CLV, aes(x = end, y = -3000, xend = start, yend = 3000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_CLV, aes(x = end, y = -2000, xend = start, yend = -2000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_CLV, aes(x = end, y = -1000, xend = start, yend = -1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_CLV, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_CLV, aes(x = end, y = 1000, xend = start, yend = 1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # Add text showing the value of each lines
  annotate("text", x = rep(max(data_CLV$id),5), y = c(-3000, -2000, -1000, 0, 1000), 
           label = c("-3000", "-2000", "-1000", "0", "1000") , 
           color="black", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=value, fill=population), 
           color = "black", size = 0.25, stat="identity", alpha=0.5) +
  ylim(-3000,1250) + theme_minimal() + coord_polar() + 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +

  # Add base line information
  geom_segment(data=base_data_CLV, aes(x = start, y = 0, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data_CLV, aes(x = title, y = 1250, label=group), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("blue", "green", "red")) +
  ggsave("FIG_3A.pdf", width = 5, height = 5)

### Creating the circle plot - Sublines without DS8
empty_bar <- 3 # Used for spacing
#### Set a number of 'empty bar' to add at the end of each group
data_sublines_noDS8 <- muts_Chrom_named_center_melt_sublines_noDS8
colnames(data_sublines_noDS8) <- c("group", "population", "value")
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data_sublines_noDS8$group), 
                             ncol(data_sublines_noDS8)) )

#### Add identifier label to dataframe
colnames(to_add) <- colnames(data_sublines_noDS8)
to_add$group <- rep(levels(data_sublines_noDS8$group), each=empty_bar)
data_sublines_noDS8 <- rbind(data_sublines_noDS8, to_add)
data_sublines_noDS8 <- data_sublines_noDS8 %>% arrange(group)
data_sublines_noDS8$id <- seq(1, nrow(data_sublines_noDS8))

# Get the name and the y position of each label
label_data_sublines_noDS8 <- data_sublines_noDS8
number_of_bar <- nrow(label_data_sublines_noDS8)
angle <- 90 - 360 * (label_data_sublines_noDS8$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data_sublines_noDS8$hjust <- ifelse( angle < -90, 1, 0)
label_data_sublines_noDS8$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a dataframe for base lines
base_data_sublines_noDS8 <- data_sublines_noDS8 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a dataframe for grid (scales)
grid_data_sublines_noDS8 <- base_data_sublines_noDS8
grid_data_sublines_noDS8$end <- grid_data_sublines_noDS8$end[ c( nrow(grid_data_sublines_noDS8), 
                                                                 1:nrow(grid_data_sublines_noDS8)-1)] + 1
grid_data_sublines_noDS8$start <- grid_data_sublines_noDS8$start - 1
grid_data_sublines_noDS8 <- grid_data_sublines_noDS8[-1,]

ggplot(data_sublines_noDS8, aes(x=as.factor(id), y=value, fill=population)) + 
  geom_bar(aes(x=as.factor(id), y=value, fill=population), stat="identity", alpha=0.5) +
  # Add reference lines
  geom_segment(data=grid_data_sublines_noDS8, aes(x = end, y = -3000, xend = start, yend = 3000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines_noDS8, aes(x = end, y = -2000, xend = start, yend = -2000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines_noDS8, aes(x = end, y = -1000, xend = start, yend = -1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines_noDS8, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines_noDS8, aes(x = end, y = 1000, xend = start, yend = 1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # Add text showing the value of each lines
  annotate("text", x = rep(max(data_sublines_noDS8$id),5), y = c(-3000, -2000, -1000, 0, 1000), 
           label = c("-3000", "-2000", "-1000", "0", "1000") , 
           color="black", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=value, fill=population), 
           color = "black", size = 0.25, stat="identity", alpha=0.5) +
  ylim(-3000,1250) + theme_minimal() + coord_polar() + 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  
  # Add base line information
  geom_segment(data=base_data_sublines_noDS8, aes(x = start, y = 0, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data_sublines_noDS8, aes(x = title, y = 1250, label=group), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "gold")) +
  ggsave("FIG_4A.pdf", width = 5, height = 5)

### Creating the circle plot - Sublines with DS8
empty_bar <- 3 # Used for spacing
#### Set a number of 'empty bar' to add at the end of each group
data_sublines <- muts_Chrom_named_center_melt_sublines
colnames(data_sublines) <- c("group", "population", "value")
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data_sublines$group), 
                             ncol(data_sublines)) )

#### Add identifier label to dataframe
colnames(to_add) <- colnames(data_sublines)
to_add$group <- rep(levels(data_sublines$group), each=empty_bar)
data_sublines <- rbind(data_sublines, to_add)
data_sublines <- data_sublines %>% arrange(group)
data_sublines$id <- seq(1, nrow(data_sublines))

# Get the name and the y position of each label
label_data_sublines <- data_sublines
number_of_bar <- nrow(label_data_sublines)
angle <- 90 - 360 * (label_data_sublines$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data_sublines$hjust <- ifelse( angle < -90, 1, 0)
label_data_sublines$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a dataframe for base lines
base_data_sublines <- data_sublines %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a dataframe for grid (scales)
grid_data_sublines <- base_data_sublines
grid_data_sublines$end <- grid_data_sublines$end[ c( nrow(grid_data_sublines), 
                                                    1:nrow(grid_data_sublines)-1)] + 1
grid_data_sublines$start <- grid_data_sublines$start - 1
grid_data_sublines <- grid_data_sublines[-1,]

ggplot(data_sublines, aes(x=as.factor(id), y=value, fill=population)) + 
  geom_bar(aes(x=as.factor(id), y=value, fill=population), stat="identity", alpha=0.5) +
  # Add reference lines
  geom_segment(data=grid_data_sublines, aes(x = end, y = -3000, xend = start, yend = 3000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines, aes(x = end, y = -2000, xend = start, yend = -2000), 
               colour = "black", alpha=1, size=0.3, inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines, aes(x = end, y = -1000, xend = start, yend = -1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data_sublines, aes(x = end, y = 1000, xend = start, yend = 1000), 
               colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # Add text showing the value of each lines
  annotate("text", x = rep(max(data_sublines$id),5), y = c(-3000, -2000, -1000, 0, 1000), 
           label = c("-3000", "-2000", "-1000", "0", "1000") , 
           color="black", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=value, fill=population), 
           color = "black", size = 0.25, stat="identity", alpha=0.5) +
  ylim(-3000,1250) + theme_minimal() + coord_polar() + 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  
  # Add base line information
  geom_segment(data=base_data_sublines, aes(x = start, y = 0, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data_sublines, aes(x = title, y = 1250, label=group), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  ggsave("FIG_S10A.pdf", width = 5, height = 5)

### COV calculation from circle plots
#### Remake dataframes by cohort for easy calculations
muts_Chrom_named_CLV <- apply(muts_Chrom_named[,c("VU", "MGH", "BR1")], 1, function(y) y)
muts_Chrom_named_melt_CLV <- melt(data = muts_Chrom_named_CLV,
                                   id.vars = rownames(muts_Chrom_named_CLV),
                                   measure.vars = colnames(muts_Chrom_named_CLV))

muts_Chrom_named_sublines_noDS8 <- apply(muts_Chrom_named[,c("DS3", "DS6", "DS7", "DS9")], 
                                         1, function(y) y)
muts_Chrom_named_melt_sublines_noDS8 <- melt(data = muts_Chrom_named_sublines_noDS8,
                                             id.vars = rownames(muts_Chrom_named_sublines_noDS8),
                                             measure.vars = colnames(muts_Chrom_named_sublines_noDS8))

muts_Chrom_named_sublines <- apply(muts_Chrom_named[,c("DS3", "DS6", "DS7", "DS8", "DS9")], 
                                         1, function(y) y)
muts_Chrom_named_melt_sublines <- melt(data = muts_Chrom_named_sublines,
                                       id.vars = rownames(muts_Chrom_named_sublines),
                                       measure.vars = colnames(muts_Chrom_named_sublines))

#### Function to calculate coefficient of variation (COV)
CV <- function(x){
  (sd(x)/mean(x))*100
}

#### COVs seen in main text figures
CLV_chromCVavg <- mean(aggregate(value ~ Var2, 
                                 data = muts_Chrom_named_melt_CLV, 
                                 FUN = CV)$value) # In FIG. 3A
sublines_noDS8_chromCVavg <- mean(aggregate(value ~ Var2, 
                                            data = muts_Chrom_named_melt_sublines_noDS8, 
                                            FUN = CV)$value) # In FIG. 4A
sublines_chromCVavg <- mean(aggregate(value ~ Var2, 
                                      data = muts_Chrom_named_melt_sublines, 
                                      FUN = CV)$value) # For reference


# Proportion and impact of unique mutations for each cohort
## Function to pull specific extra annotations into separate columns
getExtras <- function(data) {
  res_impacts <- list()
  res_classes <- list()
  res_symbols <- list()
  res_biotypes <- list()
  for (i in 1:nrow(data)) {
    res_impacts[i] <- str_match(data$Extra[i], "IMPACT=(.*?);")[,2]
    res_classes[i] <- str_match(data$Extra[i], "VARIANT_CLASS=(.*?);")[,2]
    res_symbols[i] <- str_match(data$Extra[i], "SYMBOL=(.*?);")[,2]
    res_biotypes[i] <- str_match(data$Extra[i], "BIOTYPE=(.*?)(;|$)")[,2]
  }
  data$Extra <- NULL
  data$Impact <- unlist(res_impacts)
  data$Class <- unlist(res_classes)
  data$Symbol <- unlist(res_symbols)
  data$Biotype <- unlist(res_biotypes)
  return(data)
}

# Data can be generated from VCFtoVEP.txt script or provided upon request
setwd('/Volumes/quaranta/Data/WXS/PC9/vep_data')

## Cell Line Versions
### Read in variant effect predictor (VEP) annotated variants
#### 1 = VU
#### 2 = MGH
#### 3 = BR1
sample1_vep <- read.csv('vep_sample1.txt', sep = "\t", header = T, skip = 90)
sample2_vep <- read.csv('vep_sample2.txt', sep = "\t", header = T, skip = 90)
sample3_vep <- read.csv('vep_sample3.txt', sep = "\t", header = T, skip = 90)

### Pull specific extra annotations into separate columns
test1 <- getExtras(sample1_vep)
test2 <- getExtras(sample2_vep)
test3 <- getExtras(sample3_vep)

### Remove duplicated variants (for downstream analysis)
test1_dedup <- test1[!duplicated(test1$X.Uploaded_variation), ]
test2_dedup <- test2[!duplicated(test2$X.Uploaded_variation), ]
test3_dedup <- test3[!duplicated(test3$X.Uploaded_variation), ]

### Identify intersections between cell line version variants
test_all <- Reduce(merge, list(test1, test2, test3))
test_s1_2 <- merge(test1, test2)
test_s1_2_t <- test_s1_2[which(!(unlist(test_s1_2['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]
test_s1_3 <- merge(test1, test3)
test_s1_3_t <- test_s1_3[which(!(unlist(test_s1_3['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]
test_s2_3 <- merge(test2, test3)
test_s2_3_t <- test_s2_3[which(!(unlist(test_s2_3['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]

### Identify variants unique to each cell line versions
test_s1 <- test1[which(!(unlist(test1['X.Uploaded_variation']) %in% unlist(test_s1_2_t['X.Uploaded_variation'])) &
                         !(unlist(test1['X.Uploaded_variation']) %in% unlist(test_s1_3_t['X.Uploaded_variation'])) &
                         !(unlist(test1['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]

test_s2 <- test2[which(!(unlist(test2['X.Uploaded_variation']) %in% unlist(test_s1_2_t['X.Uploaded_variation'])) &
                         !(unlist(test2['X.Uploaded_variation']) %in% unlist(test_s2_3_t['X.Uploaded_variation'])) &
                         !(unlist(test2['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]

test_s3 <- test3[which(!(unlist(test3['X.Uploaded_variation']) %in% unlist(test_s1_3_t['X.Uploaded_variation'])) &
                         !(unlist(test3['X.Uploaded_variation']) %in% unlist(test_s2_3_t['X.Uploaded_variation'])) &
                         !(unlist(test3['X.Uploaded_variation']) %in% unlist(test_all['X.Uploaded_variation']))),]

### Remove duplicated variants (variants annotations can be duplicated from VEP)
test_s1_dedup <- test_s1[!duplicated(test_s1$X.Uploaded_variation), ]
test_s2_dedup <- test_s2[!duplicated(test_s2$X.Uploaded_variation), ]
test_s3_dedup <- test_s3[!duplicated(test_s3$X.Uploaded_variation), ]

## Sublines (including DS8)
### Read in variant effect predictor (VEP) annotated variants
#### 4 = DS3
#### 5 = DS6
#### 6 = DS7
#### 7 = DS8
#### 8 = DS9
sample4_vep <- read.csv('vep_sample4.txt', sep = "\t", header = T, skip = 90)
sample5_vep <- read.csv('vep_sample5.txt', sep = "\t", header = T, skip = 90)
sample6_vep <- read.csv('vep_sample6.txt', sep = "\t", header = T, skip = 90)
sample7_vep <- read.csv('vep_sample7.txt', sep = "\t", header = T, skip = 90)
sample8_vep <- read.csv('vep_sample8.txt', sep = "\t", header = T, skip = 90)

### Pull specific extra annotations into separate columns
test4 <- getExtras(sample4_vep)
test5 <- getExtras(sample5_vep)
test6 <- getExtras(sample6_vep)
test7 <- getExtras(sample7_vep)
test8 <- getExtras(sample8_vep)

### Remove duplicated variants (for downstream analysis)
test4_dedup <- test4[!duplicated(test4$X.Uploaded_variation), ]
test5_dedup <- test5[!duplicated(test5$X.Uploaded_variation), ]
test6_dedup <- test6[!duplicated(test6$X.Uploaded_variation), ]
test7_dedup <- test7[!duplicated(test7$X.Uploaded_variation), ]
test8_dedup <- test8[!duplicated(test8$X.Uploaded_variation), ]

### Identify intersections between subline variants
test_allS <- Reduce(merge, list(test4, test5, test6, test7, test8))
test_s4_5 <- merge(test4, test5)
test_s4_6 <- merge(test4, test6)
test_s4_7 <- merge(test4, test7)
test_s4_8 <- merge(test4, test8)
test_s5_6 <- merge(test5, test6)
test_s5_7 <- merge(test5, test7)
test_s5_8 <- merge(test5, test8)
test_s6_7 <- merge(test6, test7)
test_s6_8 <- merge(test6, test8)
test_s7_8 <- merge(test7, test8)

### Identify variants unique to each subline
test_s4 <- test4[which(!(unlist(test4['X.Uploaded_variation']) %in% unlist(test_allS['X.Uploaded_variation'])) &
                       !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
                       !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
                       !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_7['X.Uploaded_variation'])) &
                       !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation']))),]
test_s5 <- test5[which(!(unlist(test5['X.Uploaded_variation']) %in% unlist(test_allS['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_7['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation']))),]
test_s6 <- test6[which(!(unlist(test6['X.Uploaded_variation']) %in% unlist(test_allS['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s6_7['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation']))),]
test_s7 <- test7[which(!(unlist(test7['X.Uploaded_variation']) %in% unlist(test_allS['X.Uploaded_variation'])) &
                         !(unlist(test7['X.Uploaded_variation']) %in% unlist(test_s4_7['X.Uploaded_variation'])) &
                         !(unlist(test7['X.Uploaded_variation']) %in% unlist(test_s5_7['X.Uploaded_variation'])) &
                         !(unlist(test7['X.Uploaded_variation']) %in% unlist(test_s6_7['X.Uploaded_variation'])) &
                         !(unlist(test7['X.Uploaded_variation']) %in% unlist(test_s7_8['X.Uploaded_variation']))),]
test_s8 <- test8[which(!(unlist(test8['X.Uploaded_variation']) %in% unlist(test_allS['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s7_8['X.Uploaded_variation']))),]

### Remove duplicated variants (variants annotations can be duplicated from VEP)
test_s4_dedup <- test_s4[!duplicated(test_s4$X.Uploaded_variation), ]
test_s5_dedup <- test_s5[!duplicated(test_s5$X.Uploaded_variation), ]
test_s6_dedup <- test_s6[!duplicated(test_s6$X.Uploaded_variation), ]
test_s7_dedup <- test_s7[!duplicated(test_s7$X.Uploaded_variation), ]
test_s8_dedup <- test_s8[!duplicated(test_s8$X.Uploaded_variation), ]


## Sublines (excluding DS8)
test_all_noDS8 <- Reduce(merge, list(test4, test5, test6, test8))
test_s4_noDS8 <- test4[which(!(unlist(test4['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
                         !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
                         !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
                         !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation']))),]
test_s5_noDS8 <- test5[which(!(unlist(test5['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
                         !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation']))),]
test_s6_noDS8 <- test6[which(!(unlist(test6['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
                         !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation']))),]
test_s8_noDS8 <- test8[which(!(unlist(test8['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation'])) &
                         !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation']))),]

### Remove duplicated variants (variants annotations can be duplicated from VEP)
test_s4_noDS8_dedup <- test_s4[!duplicated(test_s4_noDS8$X.Uploaded_variation), ]
test_s5_noDS8_dedup <- test_s5[!duplicated(test_s5_noDS8$X.Uploaded_variation), ]
test_s6_noDS8_dedup <- test_s6[!duplicated(test_s6_noDS8$X.Uploaded_variation), ]
test_s8_noDS8_dedup <- test_s8[!duplicated(test_s8_noDS8$X.Uploaded_variation), ]

## Plot unique impact mutations and proportions - CLV
### Make dataset of unique and total mutations, with the proportion
unique_prop_CLV <- data.frame(Population = c("BR1", "MGH", "VU"),
                              Total = c(nrow(test3_dedup), nrow(test2_dedup), nrow(test1_dedup)),
                              Unique = c(nrow(test_s3_dedup), nrow(test_s2_dedup), nrow(test_s1_dedup)))
unique_prop_CLV$Proportion <- unique_prop_CLV$Unique / unique_prop_CLV$Total

### Plot proportion of unique mutations in each CLV
ggplot(unique_prop_CLV, aes(x = Population, y = Proportion, fill = Population)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() + labs(x = "Cell Line Version", y = "Proportion of Unique Mutations") +
  scale_fill_manual(values = c("red", "green", "blue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(),
        plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_3B.pdf", width = 6, height = 3)

### Calculate the number of unique mutations different impact scores 
impact_VU <- dplyr::mutate(melt(table(test_s1_dedup$Impact)), prop = value / sum(value))
impact_MGH <- dplyr::mutate(melt(table(test_s2_dedup$Impact)), prop = value / sum(value))
impact_BR1 <- dplyr::mutate(melt(table(test_s3_dedup$Impact)), prop = value / sum(value))

### Put into common dataframe and remove MODIFIER variants (no predicted effect)
#### Percentage above plots in main figures are the percentage of IMPACT mutations
#### (i.e., number on this plot) in the total number of unique mutations 
impact_CLV = do.call("rbind", list(impact_VU, impact_MGH, impact_BR1))
impact_CLV$pop <- rep(c("VU", "MGH", "BR1"), each = 4)
colnames(impact_CLV) <- c("Risk", "Count", "Proportion", "Population")
impact_CLV_risky <- impact_CLV[!(impact_CLV$Risk == "MODIFIER"),]

ggplot(impact_CLV_risky, aes(x = Population, y = Count, 
                              fill = factor(Risk, levels = c("HIGH", "MODERATE", "LOW")))) + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + labs(x = "Cell Line Version", y = "Number of IMPACT Mutations") +
  scale_fill_manual(values = c("red", "orange", "yellow"), name = "Risk") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_3C.pdf", width = 6, height = 3)


## Plot unique impact mutations and proportions - sublines without DS8
### Make dataset of unique and total mutations, with the proportion
unique_prop_sublines_noDS8 <- data.frame(Population = c("DS3", "DS6", "DS7", "DS9"),
                                         Total = c(nrow(test4_dedup), nrow(test5_dedup), 
                                                   nrow(test6_dedup), nrow(test8_dedup)),
                                         Unique = c(nrow(test_s4_dedup_noDS8), nrow(test_s5_dedup_noDS8), 
                                                    nrow(test_s6_dedup_noDS8), nrow(test_s8_dedup_noDS8)))
unique_prop_sublines_noDS8$Proportion <- unique_prop_sublines_noDS8$Unique / unique_prop_sublines_noDS8$Total

### Plot proportion of unique mutations in each subline
ggplot(unique_prop_sublines_noDS8, aes(x = Population, y = Proportion, fill = Population)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() + labs(x = "Cell Line Version", y = "Proportion of Unique Mutations") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "gold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(),
        plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_4B.pdf", width = 6, height = 3)

### Calculate the number of unique mutations different impact scores 
impact_DS3_noDS8 <- dplyr::mutate(melt(table(test_s4_dedup_noDS8$Impact)), prop = value / sum(value))
impact_DS6_noDS8 <- dplyr::mutate(melt(table(test_s5_dedup_noDS8$Impact)), prop = value / sum(value))
impact_DS7_noDS8 <- dplyr::mutate(melt(table(test_s6_dedup_noDS8$Impact)), prop = value / sum(value))
impact_DS9_noDS8 <- dplyr::mutate(melt(table(test_s8_dedup_noDS8$Impact)), prop = value / sum(value))

### Put into common dataframe and remove MODIFIER variants (no predicted effect)
#### Percentage above plots in main figures are the percentage of IMPACT mutations
#### (i.e., number on this plot) in the total number of unique mutations
impact_sublines_noDS8 = do.call("rbind", list(impact_DS3_noDS8, impact_DS6_noDS8, impact_DS7_noDS8, impact_DS9_noDS8))
impact_sublines_noDS8$pop <- rep(c("DS3", "DS6", "DS7", "DS9"), each = 4)
colnames(impact_sublines_noDS8) <- c("Risk", "Count", "Proportion", "Population")
impact_sublines_noDS8_risky <- impact_sublines_noDS8[!(impact_sublines_noDS8$Risk == "MODIFIER"),]

### Plot IMPACT mutations
ggplot(impact_sublines_noDS8_risky, 
       aes(x = factor(Population, levels = c("DS3", "DS6", "DS7", "DS9")),
           y = Count, fill=factor(Risk, levels = c("HIGH", "MODERATE", "LOW")))) + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + labs(x = "Subline", y = "Number of IMPACT Mutations") +
  scale_fill_manual(values = c("red", "orange", "yellow"), name = "Risk") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12))+
  ggsave("FIG_4C.pdf", width = 6, height = 3)

### Same analysis for sublines including DS8 (highlighting DS8)
unique_prop_sublines <- data.frame(Population = c("DS3", "DS6", "DS7", "DS8", "DS9"),
                                   Total = c(nrow(test4_dedup), nrow(test5_dedup), 
                                             nrow(test6_dedup), nrow(test7_dedup),
                                             nrow(test8_dedup)),
                                   Unique = c(nrow(test_s4_dedup), nrow(test_s5_dedup), 
                                              nrow(test_s6_dedup), nrow(test_s7_dedup),
                                              nrow(test_s8_dedup)))
unique_prop_sublines$Proportion <- unique_prop_sublines$Unique / unique_prop_sublines$Total

#### Highlight DS8
unique_prop_sublines$tint <- c(0.3, 0.3, 0.3, 1, 0.3)

ggplot(unique_prop_sublines, aes(x = factor(Population, levels = c("DS3", "DS6", "DS7", "DS9", "DS8")), 
                                 y = Proportion, fill = Population, alpha = tint)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() + labs(x = "Subline", y = "Proportion of Unique Mutations") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(),
        plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_6A.pdf", width = 6, height = 3)

impact_DS3 <- dplyr::mutate(melt(table(test_s4_dedup$Impact)), prop = value / sum(value))
impact_DS6 <- dplyr::mutate(melt(table(test_s5_dedup$Impact)), prop = value / sum(value))
impact_DS7 <- dplyr::mutate(melt(table(test_s6_dedup$Impact)), prop = value / sum(value))
impact_DS8 <- dplyr::mutate(melt(table(test_s7_dedup$Impact)), prop = value / sum(value))
impact_DS9 <- dplyr::mutate(melt(table(test_s8_dedup$Impact)), prop = value / sum(value))

impact_clones = do.call("rbind", list(impact_DS3, impact_DS6, impact_DS7, impact_DS8, impact_DS9))
impact_clones$pop <- rep(c("DS3", "DS6", "DS7", "DS8", "DS9"), each = 4)
colnames(impact_clones) <- c("Risk", "Count", "Proportion", "Population")
impact_clones_risky <- impact_clones[!(impact_clones$Risk == "MODIFIER"),]

#### Highlight DS8
impact_clones_risky$tint <- rep(c(0.3, 0.3, 0.3, 1, 0.3), each = 3)

ggplot(impact_clones_risky, aes(x = factor(Population, 
                                           levels = c("DS3", "DS6", "DS7", "DS9", "DS8")),
                                y = Count, 
                                fill=factor(Risk, levels = c("HIGH", "MODERATE", "LOW")),
                                alpha = tint)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + labs(x = "Subline", y = "Number of IMPACT Mutations") +
  scale_fill_manual(values = c("red", "orange", "yellow"), name = "Risk") +
  scale_alpha_continuous(range = c(0.3, 1)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) +
  ggsave("FIG_6B.pdf", width = 6, height = 3)

# Calculate mutation class proportion
## Function to calculate the proportion of mutation classes in a pie chart
plotPie_class <- function(dat) {
  pie_dat <- table(dat)
  pie_melt <- melt(pie_dat)
  pie_melt <- dplyr::mutate(pie_melt, prop = value / sum(value))
  midpoint <- cumsum(pie_melt$prop) - pie_melt$prop/2
  pie <- ggplot(pie_melt, aes(x="", y=value, fill=factor(dat)))+
    geom_bar(stat = "identity", color = "black", width = 0.1) +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette = "BrBG") +
    guides(fill = guide_legend(title = "Group")) +
    theme_void() + theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  pie
}

## Plot the data
### CLV cohort
plotPie_class(test_s1_dedup$Class) + ggsave("FIG_3F_VU.pdf", width = 6, height = 4)
plotPie_class(test_s2_dedup$Class) + ggsave("FIG_3F_MGH.pdf", width = 6, height = 4)
plotPie_class(test_s3_dedup$Class) + ggsave("FIG_3F_BR1.pdf", width = 6, height = 4)

### Sublines (no DS8) cohort
plotPie_class(test_s4_dedup_noDS8$Class) + ggsave("FIG_4F_DS3.pdf", width = 6, height = 4)
plotPie_class(test_s5_dedup_noDS8$Class) + ggsave("FIG_4F_DS6.pdf", width = 6, height = 4)
plotPie_class(test_s6_dedup_noDS8$Class) + ggsave("FIG_4F_DS7.pdf", width = 6, height = 4)
plotPie_class(test_s8_dedup_noDS8$Class) + ggsave("FIG_4F_DS9.pdf", width = 6, height = 4)

### Sublines (with DS8) cohort
plotPie_class(test_s4_dedup$Class) + ggsave("FIG_S10D_DS3.pdf", width = 6, height = 4)
plotPie_class(test_s5_dedup$Class) + ggsave("FIG_S10D_DS6.pdf", width = 6, height = 4)
plotPie_class(test_s6_dedup$Class) + ggsave("FIG_S10D_DS7.pdf", width = 6, height = 4)
plotPie_class(test_s7_dedup$Class) + ggsave("FIG_S10D_DS8.pdf", width = 6, height = 4)
plotPie_class(test_s8_dedup$Class) + ggsave("FIG_S10D_DS9.pdf", width = 6, height = 4)


# Relatedness of SNP genotypes (hierarchical clustering, PCA)
# BiocManager::install("SNPRelate")
library(gdsfmt)
library(SNPRelate)

setwd('~/Documents/QuarantaLab/GES_2020/WES/')

## Load in file that includes all samples in common VCF file
vcf_file <- "/Users/Corey/Documents/QuarantaLab/PC9/WXS/processed_data/samples_called_vars_named.vcf"

## Convert vcf file to gds file (usable by SNPRelate)
snpgdsVCF2GDS(vcf_file,"data.gds",method ="biallelic.only")

## Open new gds file
genofile<-snpgdsOpen("data.gds")

## Run PCA analysis on SNP genotypes
library(dendextend)
set.seed(1000)

### Identify usable set of SNPs
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)

### Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

#### Calculate variance proportion (%)
pc.percent <- pca$varprop*100

#### Compile into dataframe
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#### Plot the PCA data by sample
ggplot(data = tab, aes(x = EV1, y = EV2, fill = sample.id)) + 
  geom_point(size = 4, shape = 21) + theme_bw() +
  scale_fill_manual(values = c("blue", "green", "red", "brown", 
                               "deeppink", "darkorchid", "seagreen", "gold"),
                    labels = c("VU", "MGH", "BR1", "DS3", "DS6",
                               "DS7", "DS8", "DS9"),
                    name = "Population") +
  geom_hline(yintercept=0, color = "black", size = 0.5, linetype = "dashed") +
  geom_vline(xintercept=0, color = "black", size = 0.5, linetype = "dashed") +
  labs(title = "PCA of SNP genotypes", x = "PC1", y = "PC2") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_S12A.pdf", width = 7, height = 5)

## Plot hierarchical clustering of SNP genotype relatedness
set.seed(100)

### Run HC analysis
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))

### Change labels and plot cluster tree - FIG. S12B
rv <- snpgdsCutTree(ibs.hc)
dend <- rv$dendrogram
dendextend::labels(dend) <- c("VU", "MGH", "DS8", "BR1", "DS9", "DS6", "DS3", "DS7")
plot(dend,main="Cluster Tree by SNPs - PC9", xlab = "Population", ylab = "Height")


# Identifying potential significantly mutated genes in each cohort
# install_github("im3sanger/dndscv")
library(dndscv)
library(data.table)

## Regular expression to pull out mutation characteristics
### Mutation, chromosome number, location, reference allele, mutated allele
getDatFormat <- function(data) {
  res <- str_match(data, "([^_]*)_(.*?)_(.*?)/(.*?)(/|$)")
}

## Setting the color setup for plots below
matches <- c("Essential_Splice" = "#F8766D", "Missense" = "#A3A500", "no-SNV" = "#00BF7D",
             "None" = "white", "Nonsense" = "#00B0F6", "Synonymous" = "#E76BF3")
names <- c("Essential_Splice", "Missense" , "no-SNV", "None", "Nonsense", "Synonymous")
labels <- c("Essential Splice", "Missense", "Non-SNV", "None", "Nonsense", "Synonymous")

setwd('~/Documents/QuarantaLab/GES_2020/WES/')

## Pull out mutation characteristics
muts_VU <- as.data.frame(do.call(rbind, lapply(test_s1_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_VU) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_VU$SampleID <- "VU"

muts_MGH <- as.data.frame(do.call(rbind, lapply(test_s2_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_MGH) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_MGH$SampleID <- "MGH"

muts_BR1 <- as.data.frame(do.call(rbind, lapply(test_s3_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_BR1) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_BR1$SampleID <- "BR1"

## Put mutation characteristics from cell line versions (CLV) into common dataframe
muts_CLV <- do.call("rbind", list(muts_VU, muts_MGH, muts_BR1))

## Run dNdScv (tool to identify key mutations)
### CV = NULL (because using hg38)
### Max mutations = threshold (remove hypermutator phenotype)
dndsout_CLV = dndscv(muts_CLV, refdb="RefCDS_human_GRCh38.p12.rda", 
                     cv=NULL, max_muts_per_gene_per_sample = 50)

## Pull out important mutations based on p-value threshold
CLV_keymuts <- subset(dndsout_CLV$sel_cv, pglobal_cv < 0.05) 

## Subset mutations based on important mutations
impGenes_CLV <- subset(dndsout_CLV$annotmuts, gene %in% CLV_keymuts$gene_name)

## Creating dataframe-mutation grid (for heatmap plotting)
CLV <- c("BR1", "MGH", "VU") 
fillInDF_CLV_dNdS <- expand.grid(unique(impGenes_CLV$gene), CLV)
names(fillInDF_CLV_dNdS) <- c("gene", "sampleID")
fillInDF_CLV_dNdS$interaction <- paste0(fillInDF_CLV_dNdS$gene, ":", fillInDF_CLV_dNdS$sampleID)

## Creating compiled dataframe of mutations, amino acids, and impact
CLV_all_CG_dNdS <- subset(dndsout_CLV$annotmuts, gene %in% unique(impGenes_CLV$gene))
CLV_all_CG_dNdS$interaction <- paste0(CLV_all_CG_dNdS$gene, ":", CLV_all_CG_dNdS$sampleID)

## Compile new dataframe for each mutated gene, CLV, and impact
Genes_CLV_dNdS <- merge(fillInDF_CLV_dNdS, CLV_all_CG_dNdS, c("interaction"), all.x = TRUE)
Genes_CLV_dNdS <- Genes_CLV_dNdS[,c("gene.x", "sampleID.x", "impact")]
Genes_CLV_dNdS[is.na(Genes_CLV_dNdS)] <- "None"
names(Genes_CLV_dNdS) <- c("gene", "sampleID", "impact")

## Tally number of mutations for each sample-gene pair
CG_CLV_factor_dNdS <- Genes_CLV_dNdS %>% modify_if(is.character, as.factor)
CG_CLV_factor_dNdS_count <- CG_CLV_factor_dNdS %>% 
  group_by(gene, sampleID, impact) %>% 
  tally() 

## Add shift and height column (annotations to heatmap)
CG_CLV_dNdS_factor_p <- data.table(CG_CLV_factor_dNdS_count)
CG_CLV_dNdS_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_CLV_dNdS_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Add character and numeric indecies for heatmap plotting
k <- unique(c(as.character(CG_CLV_dNdS_factor_p$sampleID)))
l <- as.numeric(factor(CG_CLV_dNdS_factor_p$sampleID, levels=k)) 

## Create gene count dataframe (for vertical heatmap bar plot)
CG_CLV_dNdS_geneCount <- subset(CG_CLV_factor_dNdS, impact != "None") %>% 
  group_by(gene, impact) %>% 
  tally() 

### Change from factor to character (for plots)
CG_CLV_dNdS_geneCount$gene <- as.character(CG_CLV_dNdS_geneCount$gene) 

## Plot vertical barplot for each gene (gene and impact classified)
ggplot(CG_CLV_dNdS_geneCount, aes(x = gene, y = n)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,20) + scale_x_discrete(limits = rev(levels(CG_CLV_dNdS_factor_p$gene))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    # axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    axis.title.x=element_blank(), axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), panel.border = element_blank()) +
  ggsave("FIG_3D_right.pdf", width = 12, height = 2)

## Plot horizontal barplot for each sample (sample and impact classified)
ggplot(impGenes_CLV %>% count(sampleID, impact), 
       aes(x = factor(sampleID, levels = c("BR1", "MGH", "VU")),
           y = n)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.125, width = 0.99) +
  theme_bw() + labs(x = "Sample", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,225) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    panel.border = element_blank()) +
  ggsave("FIG_3D_top.pdf", width = 6, height = 3)

## Add character and numeric indecies for heatmap plotting
e <- unique(c(as.character(CG_CLV_dNdS_factor_p$sampleID)))
f <- as.numeric(factor(CG_CLV_dNdS_factor_p$sampleID, levels=e))

## Plot heatmap of mutations and CLV (colored by impact and number of mutations)
### Separated bars within sample bars are if multiple mutation impacts in same CLV
ggplot(CG_CLV_dNdS_factor_p, (aes(x = gene, y = f + shift,
                                  fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Genes", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12), name = "Number of\nMutations") +
  scale_y_continuous(breaks = c(1,2,3), 
                     labels = e) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names, 
                    labels = labels) +
  theme(
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_text(size = 8),
    # axis.text.x = element_text(size = 6, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_3D_middle.pdf", width = 6, height = 12)

## Cancer-associated gene list
cancerGenes <- c("EGFR", "MET", "HER2", "BRAF", "NF1", "RAF1", "KRAS", "NRAS", "HRAS", "MAP2K2", 
                 "AKT1", "PIK3CA", "PIK3CB", "TSC1", "ALK", "APC", "EPHA2", "EPHA3", "EPHA4",
                 "ERBB4", "FGFR1", "ITK", "JAK2", "JAK3", "LRP1B", "LTK", "ROS1", "STK11", "TP53", "RB1", 
                 "ABCB5", "CFTR", "DACH1", "RELN", "CDKN2A", "DDR2", "DLEC1", "IRF1", "KEAP1", "MAP2K1",
                 "MAP3K8", "NRG1", "PPP2R1B", "PRKN", "PTEN", "RASSF1", "RET", "RIT1", "SLC22A18", "SMARCA4")

## Create new empty heatmap dataframe for cancer-associated gene mutations
fillInDF_CLV_CG <- expand.grid(cancerGenes, CLV)
names(fillInDF_CLV_CG) <- c("gene", "sampleID")
fillInDF_CLV_CG$interaction <- paste0(fillInDF_CLV_CG$gene, ":", fillInDF_CLV_CG$sampleID)

## Annotate cancer-gene associated mutation dataframe
CLV_all_CG <- subset(dndsout_CLV$annotmuts, gene %in% cancerGenes)
CLV_all_CG$interaction <- paste0(CLV_all_CG$gene, ":", CLV_all_CG$sampleID)

## Create dataframe of genes, sample name, and impact for each CLV
cancerGenes_CLV <- merge(fillInDF_CLV_CG, CLV_all_CG, c("interaction"), all.x = TRUE)
cancerGenes_CLV <- cancerGenes_CLV[,c("gene.x", "sampleID.x", "impact")]
cancerGenes_CLV[is.na(cancerGenes_CLV)] <- "None"
names(cancerGenes_CLV) <- c("gene", "sampleID", "impact")

## Convert character columns to factors
CG_CLV_factor <- cancerGenes_CLV %>% modify_if(is.character, as.factor)

## Tally the number of gene-sample-impact pairs
CG_CLV_factor <- CG_CLV_factor %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

## Calculating height and shift columns (if gene-sample pairs have more than one impact type)
CG_CLV_factor_p <- data.table(CG_CLV_factor)
CG_CLV_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_CLV_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Create character and numeric indecies (for plots below)
i <- unique(c(as.character(CG_CLV_factor_p$sampleID)))
j <- as.numeric(factor(CG_CLV_factor_p$sampleID, levels=i))

## Plot cancer-associated gene heatmap
ggplot(CG_CLV_factor_p, (aes(x = gene, y = j + shift,
                             fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Cancer Genes", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12)) +
  scale_y_continuous(breaks = c(1,2,3), 
                     labels = i) + coord_flip() +
  scale_x_discrete(limits = rev(levels(CG_CLV_factor_p$gene))) +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    # axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_3E.pdf", width = 9, height = 5)


## Generate heatmaps for sublines using same genes from CLV
### Pull out specific information from subline unique mutation lists (without DS8)
muts_DS3_noDS8 <- as.data.frame(do.call(rbind, lapply(test_s4_dedup_noDS8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS3_noDS8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS3_noDS8$SampleID <- "DS3"

muts_DS6_noDS8 <- as.data.frame(do.call(rbind, lapply(test_s5_dedup_noDS8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS6_noDS8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS6_noDS8$SampleID <- "DS6"

muts_DS7_noDS8 <- as.data.frame(do.call(rbind, lapply(test_s6_dedup_noDS8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS7_noDS8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS7_noDS8$SampleID <- "DS7"

muts_DS9_noDS8 <- as.data.frame(do.call(rbind, lapply(test_s8_dedup_noDS8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS9_noDS8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS9_noDS8$SampleID <- "DS9"

## Bind subline dataframes together
muts_sublines_noDS8 <- do.call("rbind", list(muts_DS3_noDS8, muts_DS6_noDS8, muts_DS7_noDS8, muts_DS9_noDS8))

## Run dNdScv on the subline dataframe (not to identify mutation list, but to annotate with 
## amino acids, impact, location, etc.)
dndsout_sublines_noDS8 = dndscv(muts_sublines_noDS8, refdb="RefCDS_human_GRCh38.p12.rda", 
                                cv=NULL, max_muts_per_gene_per_sample = 50)

## Pull out mutations in genes from the CLV dNdScv gene list
impGenes_sublines_noDS8_in3 <- subset(dndsout_sublines_noDS8$annotmuts, 
                                      gene %in% CLV_keymuts$gene_name)

## Generate subline list (without DS8)
sublines_noDS8 <- c("DS3", "DS6", "DS7", "DS9")

## Empty dataframe creation (same as above)
fillInDF_sublines_noDS8_in3_dNdS <- expand.grid(CLV_keymuts$gene_name, sublines_noDS8)
names(fillInDF_sublines_noDS8_in3_dNdS) <- c("gene", "sampleID")
fillInDF_sublines_noDS8_in3_dNdS$interaction <- paste0(fillInDF_sublines_noDS8_in3_dNdS$gene, ":", 
                                                       fillInDF_sublines_noDS8_in3_dNdS$sampleID)

## Subset mutations based on those in CLV list
sublines_noDS8_in3_all_CG_dNdS <- subset(dndsout_sublines_noDS8$annotmuts, 
                                         gene %in% CLV_keymuts$gene_name)
sublines_noDS8_in3_all_CG_dNdS$interaction <- paste0(sublines_noDS8_in3_all_CG_dNdS$gene, 
                                                     ":", sublines_noDS8_in3_all_CG_dNdS$sampleID)

## Create compiled dataframe of gene, sample, and impact
genes_sublines_noDS8_in3_dNdS <- merge(fillInDF_sublines_noDS8_in3_dNdS, 
                                       sublines_noDS8_in3_all_CG_dNdS, c("interaction"), all.x = TRUE)
genes_sublines_noDS8_in3_dNdS <- genes_sublines_noDS8_in3_dNdS[,c("gene.x", "sampleID.x", "impact")]
genes_sublines_noDS8_in3_dNdS[is.na(genes_sublines_noDS8_in3_dNdS)] <- "None"
names(genes_sublines_noDS8_in3_dNdS) <- c("gene", "sampleID", "impact")

## Convert character columns to factors
CG_sublines_noDS8_in3_factor_dNdS <- genes_sublines_noDS8_in3_dNdS %>% modify_if(is.character, as.factor)

## Tally gene-sample-impact pairs
CG_sublines_noDS8_in3_factor_dNdS_count <- CG_sublines_noDS8_in3_factor_dNdS %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

## Creating columns to show multiple mutation impacts in each sample-gene pair
CG_sublines_noDS8_in3_dNdS_factor_p <- data.table(CG_sublines_noDS8_in3_factor_dNdS_count)
CG_sublines_noDS8_in3_dNdS_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_sublines_noDS8_in3_dNdS_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Modifying dataframe to assist in counts (below)
CG_sublines_noDS8_in3_dNdS_factor_p$gene <- as.character(CG_sublines_noDS8_in3_dNdS_factor_p$gene)
CG_sublines_noDS8_dNdScvCLV_geneCount <- within(CG_sublines_noDS8_in3_dNdS_factor_p, n[n == 1 & impact == "None"] <- 0)
CG_sublines_noDS8_dNdScvCLV_geneCount <- CG_sublines_noDS8_dNdScvCLV_geneCount[, .(count = .N, var = sum(n)), by = list(gene, impact)]
CG_sublines_noDS8_dNdScvCLV_geneCount <- within(CG_sublines_noDS8_dNdScvCLV_geneCount, count[impact == "None"] <- 0)

## Plot gene count barplot (vertical) for sublines (annotated by impact)
ggplot(CG_sublines_noDS8_dNdScvCLV_geneCount, aes(x = gene, y = count)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,20) + scale_x_discrete(limits = rev(levels(CG_CLV_dNdS_factor_p$gene))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    axis.title.x=element_blank(), axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), panel.border = element_blank()) +
  ggsave("FIG_4D_right.pdf", width = 12, height = 2)

## Plot subline count barplot (horizontal) for genes (annotated by impact)
ggplot(impGenes_sublines_noDS8_in3 %>% count(sampleID, impact), 
       aes(x = factor(sampleID, levels = c("DS3", "DS6", "DS7", "DS9")),
           y = n)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.125, width = 0.99) +
  theme_bw() + labs(x = "Sample", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,225) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    panel.border = element_blank()) +
  ggsave("FIG_4D_top.pdf", width = 6, height = 3)

## Create character and numeric indecies (for below plot)
a <- unique(c(as.character(CG_sublines_noDS8_in3_dNdS_factor_p$sampleID)))
b <- as.numeric(factor(CG_sublines_noDS8_in3_dNdS_factor_p$sampleID, levels=a))

## Plot subline mutations for CLV distinguishing genes (colored by impact and number of mutations)
ggplot(CG_sublines_noDS8_in3_dNdS_factor_p, aes(x = gene, y = b + shift,
                                                fill = impact, alpha = n, height = height)) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Genes (Distinguishing Cell Line Versions)", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12), name = "Number of\nMutations") +
  scale_y_continuous(breaks = c(1,2,3,4), 
                     labels = a) + 
  scale_x_discrete(limits = levels(CG_CLV_dNdS_factor_p$gene)) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_4D_middle.pdf", width = 6, height = 12)

## Create cancer-associated gene signature dataframe for sublines (without DS8)
fillInDF_sublines_noDS8_in3 <- expand.grid(cancerGenes, sublines_noDS8)
names(fillInDF_sublines_noDS8_in3) <- c("gene", "sampleID")
fillInDF_sublines_noDS8_in3$interaction <- paste0(fillInDF_sublines_noDS8_in3$gene, ":", 
                                                  fillInDF_sublines_noDS8_in3$sampleID)

## Annotate cancer-associated mutation data by gene and sample (for below)
sublines_noDS8_in3_all_CG <- subset(dndsout_sublines_noDS8$annotmuts, gene %in% cancerGenes)
sublines_noDS8_in3_all_CG$interaction <- paste0(sublines_noDS8_in3_all_CG$gene, 
                                                ":", sublines_noDS8_in3_all_CG$sampleID)

## Annotate signature dataframe with mutation data
cancerGenes_sublines_noDS8_in3 <- merge(fillInDF_sublines_noDS8_in3, sublines_noDS8_in3_all_CG, 
                                        c("interaction"), all.x = TRUE)
cancerGenes_sublines_noDS8_in3 <- cancerGenes_sublines_noDS8_in3[,c("gene.x", "sampleID.x", "impact")]
cancerGenes_sublines_noDS8_in3[is.na(cancerGenes_sublines_noDS8_in3)] <- "None"
names(cancerGenes_sublines_noDS8_in3) <- c("gene", "sampleID", "impact")

## Modify character columns to factors (for below)
CG_sublines_noDS8_in3_factor <- cancerGenes_sublines_noDS8_in3 %>% modify_if(is.character, as.factor)

## Count number of gene-sample-impact pairs
CG_sublines_noDS8_in3_factor <- CG_sublines_noDS8_in3_factor %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

## Create column for spacing in heatmap (if multiple impacts for each gene-sample pair)
CG_sublines_noDS8_in3_factor_p <- data.table(CG_sublines_noDS8_in3_factor)
CG_sublines_noDS8_in3_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_sublines_noDS8_in3_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Add indecies (for plotting)
c <- unique(c(as.character(CG_sublines_noDS8_in3_factor_p$sampleID)))
d <- as.numeric(factor(CG_sublines_noDS8_in3_factor_p$sampleID, levels=c))

ggplot(CG_sublines_noDS8_in3_factor_p, (aes(x = gene, y = d + shift,
                                            fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(CG_CLV_factor_p$gene))) +
  labs(x = "Cancer Genes", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12)) +
  scale_y_continuous(breaks = c(1,2,3,4), 
                     labels = c) +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_4E.pdf", width = 9, height = 5)


## Create mutation lists of sublines (including DS8, for input into dNdScv)
muts_DS3 <- as.data.frame(do.call(rbind, lapply(test_s4_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS3) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS3$SampleID <- "DS3"

muts_DS6 <- as.data.frame(do.call(rbind, lapply(test_s5_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS6) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS6$SampleID <- "DS6"

muts_DS7 <- as.data.frame(do.call(rbind, lapply(test_s6_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS7) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS7$SampleID <- "DS7"

muts_DS8 <- as.data.frame(do.call(rbind, lapply(test_s7_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS8$SampleID <- "DS8"

muts_DS9 <- as.data.frame(do.call(rbind, lapply(test_s8_dedup$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS9) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS9$SampleID <- "DS9"

## Bind sublines together
muts_sublines <- do.call("rbind", list(muts_DS3, muts_DS6, muts_DS7, muts_DS8, muts_DS9))

## Run dNdScv on sublines (just to annotate mutations)
dndsout_sublines = dndscv(muts_sublines, refdb="RefCDS_human_GRCh38.p12.rda", 
                          cv=NULL, max_muts_per_gene_per_sample = 50)

## Subset mutations by CLV mutation list
impGenes_sublines_in3 <- subset(dndsout_sublines$annotmuts, gene %in% CLV_keymuts$gene_name)

## Enumerate subline list
sublines <- c("DS3", "DS6", "DS7", "DS8", "DS9")

## Create empty dataframe 
fillInDF_sublines_in3_dNdS <- expand.grid(CLV_keymuts$gene_name, sublines)
names(fillInDF_sublines_in3_dNdS) <- c("gene", "sampleID")
fillInDF_sublines_in3_dNdS$interaction <- paste0(fillInDF_sublines_in3_dNdS$gene, ":", 
                                                 fillInDF_sublines_in3_dNdS$sampleID)

## Annotate mutation data with sample-gene pair
sublines_in3_all_CG_dNdS <- subset(dndsout_sublines$annotmuts, gene %in% CLV_keymuts$gene_name)
sublines_in3_all_CG_dNdS$interaction <- paste0(sublines_in3_all_CG_dNdS$gene, ":", sublines_in3_all_CG_dNdS$sampleID)

## Merge two datasframes to generate new dataframe
genes_sublines_in3_dNdS <- merge(fillInDF_sublines_in3_dNdS, 
                                 sublines_in3_all_CG_dNdS, c("interaction"), all.x = TRUE)
genes_sublines_in3_dNdS <- genes_sublines_in3_dNdS[,c("gene.x", "sampleID.x", "impact")]
genes_sublines_in3_dNdS[is.na(genes_sublines_in3_dNdS)] <- "None"
names(genes_sublines_in3_dNdS) <- c("gene", "sampleID", "impact")

## Change character columns into factors
CG_sublines_in3_factor_dNdS <- genes_sublines_in3_dNdS %>% modify_if(is.character, as.factor)

## Count the number of gene-sample-impact pairs
CG_sublines_in3_factor_dNdS_count <- CG_sublines_in3_factor_dNdS %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

## Add columns that split heatmap cells if multiple impact types
CG_sublines_in3_dNdS_factor_p <- data.table(CG_sublines_in3_factor_dNdS_count)
CG_sublines_in3_dNdS_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_sublines_in3_dNdS_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Modify dataframe so it's easier to plot below
CG_sublines_in3_dNdS_factor_p$gene <- as.character(CG_sublines_in3_dNdS_factor_p$gene)
CG_sublines_dNdScvCLV_geneCount <- within(CG_sublines_in3_dNdS_factor_p, n[n == 1 & impact == "None"] <- 0)
CG_sublines_dNdScvCLV_geneCount <- CG_sublines_dNdScvCLV_geneCount[, .(count = .N, var = sum(n)), by = list(gene, impact)]
CG_sublines_dNdScvCLV_geneCount <- within(CG_sublines_dNdScvCLV_geneCount, count[impact == "None"] <- 0)

## Plot mutation count for each gene for sublines (inculding DS8)
ggplot(CG_sublines_dNdScvCLV_geneCount, aes(x = gene, y = count)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,20) + scale_x_discrete(limits = rev(levels(CG_CLV_dNdS_factor_p$gene))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    # axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    axis.title.x=element_blank(), axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), panel.border = element_blank()) +
  ggsave("FIG_S10B_right.pdf", width = 12, height = 2)

## Plot mutation count for each subline across all genes (including DS8)
ggplot(impGenes_sublines_in3 %>% count(sampleID, impact), 
       aes(x = factor(sampleID, levels = c("DS3", "DS6", "DS7", "DS8", "DS9")),
           y = n)) +
  geom_bar(stat = "identity", color = "black", aes(fill = impact), size = 0.125, width = 0.99) +
  theme_bw() + labs(x = "Sample", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  ylim(0,225) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=8),
    legend.title = element_text(size=12), axis.title=element_text(size=10),
    panel.border = element_blank()) +
  ggsave("FIG_S10B_top.pdf", width = 6, height = 3)

## Identify indecies for plots below
w <- unique(c(as.character(CG_sublines_in3_dNdS_factor_p$sampleID)))
x <- as.numeric(factor(CG_sublines_in3_dNdS_factor_p$sampleID, levels=w))

## Plot subline (including DS8) heatmap (CLV genes)
ggplot(CG_sublines_in3_dNdS_factor_p, (aes(x = gene, y = x + shift,
                                           fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Genes (Distinguishing Cell Line Versions)", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12), name = "Number of\nMutations") +
  scale_y_continuous(breaks = c(1,2,3,4,5), 
                     labels = w) + 
  scale_x_discrete(limits = levels(CG_CLV_dNdS_factor_p$gene)) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(), 
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    # axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_S10B_middle.pdf", width = 6, height = 12)

## Calculate cancer-associated heatmap for sublines (including DS8)
### Create open dataframe
fillInDF_sublines_in3 <- expand.grid(cancerGenes, sublines)
names(fillInDF_sublines_in3) <- c("gene", "sampleID")
fillInDF_sublines_in3$interaction <- paste0(fillInDF_sublines_in3$gene, 
                                            ":", fillInDF_sublines_in3$sampleID)

## Subset for only cancer genes
sublines_in3_all_CG <- subset(dndsout_sublines$annotmuts, gene %in% cancerGenes)
sublines_in3_all_CG$interaction <- paste0(sublines_in3_all_CG$gene, ":", sublines_in3_all_CG$sampleID)

## Annotate open dataframe with cancer gene information
cancerGenes_sublines_in3 <- merge(fillInDF_sublines_in3, sublines_in3_all_CG, c("interaction"), all.x = TRUE)
cancerGenes_sublines_in3 <- cancerGenes_sublines_in3[,c("gene.x", "sampleID.x", "impact")]
cancerGenes_sublines_in3[is.na(cancerGenes_sublines_in3)] <- "None"
# cancerGenes_clones_noDS8$n[is.na(cancerGenes_clones_noDS8$n)] <- 0
# cancerGenes_clones_noDS8 <- cancerGenes_clones_noDS8[,c(2,3,6)]
names(cancerGenes_sublines_in3) <- c("gene", "sampleID", "impact")

## Convert character columns to factors
CG_sublines_in3_factor <- cancerGenes_sublines_in3 %>% modify_if(is.character, as.factor)

## Count the number of gene-sample-impact pairs
CG_sublines_in3_factor <- CG_sublines_in3_factor %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

## Add shift and height columns (for positioning on heatmap below)
CG_sublines_in3_factor_p <- data.table(CG_sublines_in3_factor)
CG_sublines_in3_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_sublines_in3_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Create indecies for plots below
y <- unique(c(as.character(CG_sublines_in3_factor_p$sampleID)))
z <- as.numeric(factor(CG_sublines_in3_factor_p$sampleID, levels=y))

ggplot(CG_sublines_in3_factor_p, (aes(x = gene, y = z + shift,
                                      fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Cancer Genes", y = "Sample") +
  scale_alpha(range = c(0.5,1), limits = c(1,12)) +
  scale_y_continuous(breaks = c(1,2,3,4,5), 
                     labels = y) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(CG_CLV_factor_p$gene))) +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12)) + 
  ggsave("FIG_S10C.pdf", width = 5, height = 7)

# Save data for GO semantic similarity analysis
save(test_s1_dedup, test_s2_dedup, test_s3_dedup,
     test_s4_dedup_noDS8, test_s5_dedup_noDS8, test_s6_dedup_noDS8, test_s8_dedup_noDS8,
     test_s4_dedup, test_s5_dedup, test_s6_dedup, test_s7_dedup, test_s8_dedup,
     file = "../GO/variants_byCohort.RData")