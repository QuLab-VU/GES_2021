setwd('~/git/GES_2020/WES/')

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
  ggsave("FIG_S3A.pdf", width = 12, height = 4)

# Initiate objects for data quality control (QC) analysis 
vcf <- read.vcfR("~/git/GES_2020/WES/samples_called_vars_named.vcf.gz", verbose = TRUE)

# Download from: https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large
dna <- ape::read.dna("Homo_sapiens_assembly38.fasta", format = "fasta")

# Download from: ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/
## I pulled directly from removed genes.gtf file that is the same as above
## File used to be located here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
gff <- read.table("genes.gtf", sep="\t", quote="")

# Plot QC metrics for compiled dataset (all variants in all samples) - FIG. S3B
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

### FIG. S3C - saved as 6x9 PDF
upset(fromExpression(upset_CLV), order.by = "freq",
      sets.bar.color = c("green", "red", "blue"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 101000,
      mainbar.y.label = "Number of Mutations", sets.x.label = "Set Size", set_size.scale_max = 160000,
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

### Sublines
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

### FIG. S3D - saved as 6x9 PDF
upset(fromExpression(upset_sublines), order.by = "freq",
      sets.bar.color = c("brown", "seagreen", "deeppink", "gold", "darkorchid"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 106000,
      mainbar.y.label = "Number of Mutations", sets.x.label = "Set Size", set_size.scale_max = 160000,
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))


### Sublines including VU
shared_vars_VUDS <- read.csv("shared_variants_VUDSlines.csv", header = T)
upset_VUDS <- c("VU" = shared_vars_VUDS[1, "Number"],
                "DS3" = shared_vars_VUDS[2, "Number"],
                "DS6" = shared_vars_VUDS[3, "Number"],
                "DS7" = shared_vars_VUDS[4, "Number"],
                "DS8" = shared_vars_VUDS[5, "Number"],
                "DS9" = shared_vars_VUDS[6, "Number"],
                "VU&DS3" = shared_vars_VUDS[7, "Number"],
                "VU&DS6" = shared_vars_VUDS[8, "Number"],
                "VU&DS7" = shared_vars_VUDS[9, "Number"],
                "VU&DS8" = shared_vars_VUDS[10, "Number"],
                "VU&DS9" = shared_vars_VUDS[11, "Number"],
                "DS3&DS6" = shared_vars_VUDS[12, "Number"],
                "DS3&DS7" = shared_vars_VUDS[13, "Number"],
                "DS3&DS8" = shared_vars_VUDS[14, "Number"],
                "DS3&DS9" = shared_vars_VUDS[15, "Number"],
                "DS6&DS7" = shared_vars_VUDS[16, "Number"],
                "DS6&DS8" = shared_vars_VUDS[17, "Number"],
                "DS6&DS9" = shared_vars_VUDS[18, "Number"],
                "DS7&DS8" = shared_vars_VUDS[19, "Number"],
                "DS7&DS9" = shared_vars_VUDS[20, "Number"],
                "DS8&DS9" = shared_vars_VUDS[21, "Number"],
                "VU&DS3&DS6" = shared_vars_VUDS[22, "Number"],
                "VU&DS3&DS7" = shared_vars_VUDS[23, "Number"],
                "VU&DS3&DS8" = shared_vars_VUDS[24, "Number"],
                "VU&DS3&DS9" = shared_vars_VUDS[25, "Number"],
                "VU&DS6&DS7" = shared_vars_VUDS[26, "Number"],
                "VU&DS6&DS8" = shared_vars_VUDS[27, "Number"],
                "VU&DS6&DS9" = shared_vars_VUDS[28, "Number"],
                "VU&DS7&DS8" = shared_vars_VUDS[29, "Number"],
                "VU&DS7&DS9" = shared_vars_VUDS[30, "Number"],
                "VU&DS8&DS9" = shared_vars_VUDS[31, "Number"],
                "DS3&DS6&DS7" = shared_vars_VUDS[32, "Number"],
                "DS3&DS6&DS8" = shared_vars_VUDS[33, "Number"],
                "DS3&DS6&DS9" = shared_vars_VUDS[34, "Number"],
                "DS3&DS7&DS8" = shared_vars_VUDS[35, "Number"],
                "DS3&DS7&DS9" = shared_vars_VUDS[36, "Number"],
                "DS3&DS8&DS9" = shared_vars_VUDS[37, "Number"],
                "DS6&DS7&DS8" = shared_vars_VUDS[38, "Number"],
                "DS6&DS7&DS9" = shared_vars_VUDS[39, "Number"],
                "DS6&DS8&DS9" = shared_vars_VUDS[40, "Number"],
                "DS7&DS8&DS9" = shared_vars_VUDS[41, "Number"],
                "VU&DS3&DS6&DS7" = shared_vars_VUDS[42, "Number"],
                "VU&DS3&DS6&DS8" = shared_vars_VUDS[43, "Number"],
                "VU&DS3&DS6&DS9" = shared_vars_VUDS[44, "Number"],
                "VU&DS3&DS7&DS8" = shared_vars_VUDS[45, "Number"],
                "VU&DS3&DS7&DS9" = shared_vars_VUDS[46, "Number"],
                "VU&DS3&DS8&DS9" = shared_vars_VUDS[47, "Number"],
                "VU&DS6&DS7&DS8" = shared_vars_VUDS[48, "Number"],
                "VU&DS6&DS7&DS9" = shared_vars_VUDS[49, "Number"],
                "VU&DS6&DS8&DS9" = shared_vars_VUDS[50, "Number"],
                "VU&DS7&DS8&DS9" = shared_vars_VUDS[51, "Number"],
                "DS3&DS6&DS7&DS8" = shared_vars_VUDS[52, "Number"],
                "DS3&DS6&DS7&DS9" = shared_vars_VUDS[53, "Number"],
                "DS3&DS6&DS8&DS9" = shared_vars_VUDS[54, "Number"],
                "DS3&DS7&DS8&DS9" = shared_vars_VUDS[55, "Number"],
                "DS6&DS7&DS8&DS9" = shared_vars_VUDS[56, "Number"],
                "VU&DS3&DS6&DS7&DS8" = shared_vars_VUDS[57, "Number"],
                "VU&DS3&DS6&DS7&DS9" = shared_vars_VUDS[58, "Number"],
                "VU&DS3&DS6&DS8&DS9" = shared_vars_VUDS[59, "Number"],
                "VU&DS3&DS7&DS8&DS9" = shared_vars_VUDS[60, "Number"],
                "VU&DS6&DS7&DS8&DS9" = shared_vars_VUDS[61, "Number"],
                "DS3&DS6&DS7&DS8&DS9" = shared_vars_VUDS[62, "Number"],
                "VU&DS3&DS6&DS7&DS8&DS9" = shared_vars_VUDS[63, "Number"])

### FIG. S3E - saved as 6x9 PDF
upset(fromExpression(upset_VUDS), order.by = "freq", nsets = 6,
      sets.bar.color = c("brown", "seagreen", "deeppink", "gold", "darkorchid", "blue"),
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
#### Sublines 
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

ggplot() + 
  geom_bar(data = data_CLV, aes(x=as.factor(id), y=value, fill=population), color = "black",
           stat="identity", alpha=0.5) +
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
           label = c("-3000", "-2000", "-1000", "0", "1000"), 
           color="black", size=3 , angle=0, fontface="bold", hjust=1) +
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

### Creating the circle plot - Sublines 
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

ggplot() + 
  geom_bar(data_sublines, aes(x=as.factor(id), y=value, fill=population), 
           color = "black", stat="identity", alpha=0.5) +
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
  ggsave("FIG_4A_REVISED.pdf", width = 5, height = 5)

### COV calculation from circle plots
#### Remake dataframes by cohort for easy calculations
muts_Chrom_named_CLV <- apply(muts_Chrom_named[,c("VU", "MGH", "BR1")], 1, function(y) y)
muts_Chrom_named_melt_CLV <- melt(data = muts_Chrom_named_CLV,
                                   id.vars = rownames(muts_Chrom_named_CLV),
                                   measure.vars = colnames(muts_Chrom_named_CLV))

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
sublines_chromCVavg <- mean(aggregate(value ~ Var2, 
                                      data = muts_Chrom_named_melt_sublines, 
                                      FUN = CV)$value) # In revised FIG. 4A


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

# Data can be generated from VCFtoVEP.txt script or in supplementary data files
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


# ## Sublines (excluding DS8)
# test_all_noDS8 <- Reduce(merge, list(test4, test5, test6, test8))
# test_s4_noDS8 <- test4[which(!(unlist(test4['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
#                          !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
#                          !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
#                          !(unlist(test4['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation']))),]
# test_s5_noDS8 <- test5[which(!(unlist(test5['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
#                          !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s4_5['X.Uploaded_variation'])) &
#                          !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
#                          !(unlist(test5['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation']))),]
# test_s6_noDS8 <- test6[which(!(unlist(test6['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
#                          !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s4_6['X.Uploaded_variation'])) &
#                          !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s5_6['X.Uploaded_variation'])) &
#                          !(unlist(test6['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation']))),]
# test_s8_noDS8 <- test8[which(!(unlist(test8['X.Uploaded_variation']) %in% unlist(test_all_noDS8['X.Uploaded_variation'])) &
#                          !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s4_8['X.Uploaded_variation'])) &
#                          !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s5_8['X.Uploaded_variation'])) &
#                          !(unlist(test8['X.Uploaded_variation']) %in% unlist(test_s6_8['X.Uploaded_variation']))),]
# 
# ### Remove duplicated variants (variants annotations can be duplicated from VEP)
# test_s4_noDS8_dedup <- test_s4[!duplicated(test_s4_noDS8$X.Uploaded_variation), ]
# test_s5_noDS8_dedup <- test_s5[!duplicated(test_s5_noDS8$X.Uploaded_variation), ]
# test_s6_noDS8_dedup <- test_s6[!duplicated(test_s6_noDS8$X.Uploaded_variation), ]
# test_s8_noDS8_dedup <- test_s8[!duplicated(test_s8_noDS8$X.Uploaded_variation), ]

## Plot unique impact mutations and proportions - CLV
### Make dataset of unique and total mutations, with the proportion
unique_prop_CLV <- data.frame(Population = c("BR1", "MGH", "VU"),
                              Total = c(nrow(test3_dedup), nrow(test2_dedup), nrow(test1_dedup)),
                              Unique = c(nrow(test_s3_dedup), nrow(test_s2_dedup), nrow(test_s1_dedup)))
unique_prop_CLV$Proportion <- unique_prop_CLV$Unique / unique_prop_CLV$Total

### Plot proportion of unique mutations in each CLV
unique_prop_CLV$Population <- factor(unique_prop_CLV$Population,
                                     levels = c("BR1", "MGH", "VU"))
ggplot(unique_prop_CLV, aes(x = Population, y = Proportion, fill = Population)) + 
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() + labs(x = "Cell Line Version", y = "Proportion of Unique Mutations") +
  scale_fill_manual(values = c("VU" = "blue", "MGH" = "green", "BR1" = "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(), axis.title = element_text(size = 14),
        axis.text = element_text(size = 14), legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
ggsave("FIG_3B_REVISION.pdf", width = 6, height = 3)

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

impact_CLV_risky$Population <- factor(impact_CLV_risky$Population,
                                      levels = c("BR1", "MGH", "VU"))
impact_CLV_risky$Risk <- factor(impact_CLV_risky$Risk,
                                levels = c("HIGH", "MODERATE", "LOW"))
ggplot(impact_CLV_risky, aes(x = Population, y = Count, fill = Risk)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + labs(x = "Subline", y = "Number of IMPACT Mutations") +
  scale_fill_manual(values = c("HIGH" = "red", "MODERATE" = "orange", "LOW" = "yellow"), 
                    name = "Risk") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.text=element_text(size=14),
    axis.title=element_text(size=14)) +
  ggsave("FIG_3C_REVISION.pdf", width = 6, height = 3)


# ## Plot unique impact mutations and proportions - sublines without DS8
# ### Make dataset of unique and total mutations, with the proportion
# unique_prop_sublines_noDS8 <- data.frame(Population = c("DS3", "DS6", "DS7", "DS9"),
#                                          Total = c(nrow(test4_dedup), nrow(test5_dedup), 
#                                                    nrow(test6_dedup), nrow(test8_dedup)),
#                                          Unique = c(nrow(test_s4_dedup_noDS8), nrow(test_s5_dedup_noDS8), 
#                                                     nrow(test_s6_dedup_noDS8), nrow(test_s8_dedup_noDS8)))
# unique_prop_sublines_noDS8$Proportion <- unique_prop_sublines_noDS8$Unique / unique_prop_sublines_noDS8$Total
# 
# ### Plot proportion of unique mutations in each subline
# ggplot(unique_prop_sublines_noDS8, aes(x = Population, y = Proportion, fill = Population)) +
#   geom_bar(stat = "identity", color = "black") + 
#   theme_minimal() + labs(x = "Cell Line Version", y = "Proportion of Unique Mutations") +
#   scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "gold")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(), axis.line.y = element_line(),
#         plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#         legend.title = element_text(size=12), axis.title=element_text(size=12)) +
#   ggsave("FIG_4B.pdf", width = 6, height = 3)
# 
# ### Calculate the number of unique mutations different impact scores 
# impact_DS3_noDS8 <- dplyr::mutate(melt(table(test_s4_dedup_noDS8$Impact)), prop = value / sum(value))
# impact_DS6_noDS8 <- dplyr::mutate(melt(table(test_s5_dedup_noDS8$Impact)), prop = value / sum(value))
# impact_DS7_noDS8 <- dplyr::mutate(melt(table(test_s6_dedup_noDS8$Impact)), prop = value / sum(value))
# impact_DS9_noDS8 <- dplyr::mutate(melt(table(test_s8_dedup_noDS8$Impact)), prop = value / sum(value))
# 
# ### Put into common dataframe and remove MODIFIER variants (no predicted effect)
# #### Percentage above plots in main figures are the percentage of IMPACT mutations
# #### (i.e., number on this plot) in the total number of unique mutations
# impact_sublines_noDS8 = do.call("rbind", list(impact_DS3_noDS8, impact_DS6_noDS8, impact_DS7_noDS8, impact_DS9_noDS8))
# impact_sublines_noDS8$pop <- rep(c("DS3", "DS6", "DS7", "DS9"), each = 4)
# colnames(impact_sublines_noDS8) <- c("Risk", "Count", "Proportion", "Population")
# impact_sublines_noDS8_risky <- impact_sublines_noDS8[!(impact_sublines_noDS8$Risk == "MODIFIER"),]
# 
# ### Plot IMPACT mutations
# ggplot(impact_sublines_noDS8_risky, 
#        aes(x = factor(Population, levels = c("DS3", "DS6", "DS7", "DS9")),
#            y = Count, fill=factor(Risk, levels = c("HIGH", "MODERATE", "LOW")))) + 
#   geom_bar(stat = "identity", color = "black") + 
#   theme_classic() + labs(x = "Subline", y = "Number of IMPACT Mutations") +
#   scale_fill_manual(values = c("red", "orange", "yellow"), name = "Risk") +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#     legend.title = element_text(size=12), axis.title=element_text(size=12))+
#   ggsave("FIG_4C.pdf", width = 6, height = 3)

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
# unique_prop_sublines$tint <- c(0.3, 0.3, 0.3, 1, 0.3)

ggplot(unique_prop_sublines, aes(x = factor(Population, levels = c("DS3", "DS6", "DS7", "DS9", "DS8")), 
                                 y = Proportion, fill = Population)) + 
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() + labs(x = "Subline", y = "Proportion of Unique Mutations") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  # scale_alpha_continuous(range = c(0.3, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(), axis.title = element_text(size = 14),
        axis.text = element_text(size = 14), legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggsave("FIG_4B_REVISION.pdf", width = 6, height = 3)


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
# impact_clones_risky$tint <- rep(c(0.3, 0.3, 0.3, 1, 0.3), each = 3)

ggplot(impact_clones_risky, aes(x = factor(Population, 
                                           levels = c("DS3", "DS6", "DS7", "DS9", "DS8")),
                                y = Count, 
                                fill=factor(Risk, levels = c("HIGH", "MODERATE", "LOW")))) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + labs(x = "Subline", y = "Number of IMPACT Mutations") +
  scale_fill_manual(values = c("red", "orange", "yellow"), name = "Risk") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.text=element_text(size=14),
    axis.title=element_text(size=14)) +
  ggsave("FIG_4C_REVISION.pdf", width = 6, height = 3)

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
plotPie_class(test_s1_dedup$Class) + ggsave("FIG_S4B_VU.pdf", width = 6, height = 4)
plotPie_class(test_s2_dedup$Class) + ggsave("FIG_S4B_MGH.pdf", width = 6, height = 4)
plotPie_class(test_s3_dedup$Class) + ggsave("FIG_S4B_BR1.pdf", width = 6, height = 4)

# ### Sublines (no DS8) cohort
# plotPie_class(test_s4_dedup_noDS8$Class) + ggsave("FIG_4F_DS3.pdf", width = 6, height = 4)
# plotPie_class(test_s5_dedup_noDS8$Class) + ggsave("FIG_4F_DS6.pdf", width = 6, height = 4)
# plotPie_class(test_s6_dedup_noDS8$Class) + ggsave("FIG_4F_DS7.pdf", width = 6, height = 4)
# plotPie_class(test_s8_dedup_noDS8$Class) + ggsave("FIG_4F_DS9.pdf", width = 6, height = 4)

### Sublines (with DS8) cohort
plotPie_class(test_s4_dedup$Class) + ggsave("FIG_S4B_DS3.pdf", width = 6, height = 4)
plotPie_class(test_s5_dedup$Class) + ggsave("FIG_S4B_DS6.pdf", width = 6, height = 4)
plotPie_class(test_s6_dedup$Class) + ggsave("FIG_S4B_DS7.pdf", width = 6, height = 4)
plotPie_class(test_s7_dedup$Class) + ggsave("FIG_S4B_DS8.pdf", width = 6, height = 4)
plotPie_class(test_s8_dedup$Class) + ggsave("FIG_S4B_DS9.pdf", width = 6, height = 4)


# Relatedness of SNP genotypes (hierarchical clustering, PCA)
# BiocManager::install("SNPRelate")
library(gdsfmt)
library(SNPRelate)

setwd('~/git/GES_2020/WES/')

## Load in file that includes all samples in common VCF file - see supplementary information files
vcf_file <- "/Users/Corey/Documents/QuarantaLab/PC9/WXS/accre_vars/samples_called_vars_named.vcf"

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
  ggsave("FIG_S14A.pdf", width = 7, height = 5)

## Plot hierarchical clustering of SNP genotype relatedness
set.seed(100)

### Run HC analysis
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))

### Change labels and plot cluster tree - FIG. S14B
rv <- snpgdsCutTree(ibs.hc)
dend <- rv$dendrogram
dendextend::labels(dend) <- c("VU", "MGH", "DS8", "BR1", "DS9", "DS6", "DS3", "DS7")
plot(dend,main="Cluster Tree by SNPs - PC9", xlab = "Population", ylab = "Height")


# Identifying potential significantly mutated genes in each cohort
## Find subsets unique among all 8 cell populations, instead of CLV and subline cohorts
test1_dedup$Population <- NULL
test2_dedup$Population <- NULL
test3_dedup$Population <- NULL
test4_dedup$Population <- NULL
test5_dedup$Population <- NULL
test6_dedup$Population <- NULL
test7_dedup$Population <- NULL
test8_dedup$Population <- NULL

test_all_8 <- Reduce(merge, list(test1_dedup, test2_dedup, test3_dedup,
                                 test4_dedup, test5_dedup, test6_dedup,
                                 test7_dedup, test8_dedup))
# 2 way comparisons - 28 instances - used for removing non-unique mutations
test_dedup_s1_2 <- merge(test1_dedup, test2_dedup)
test_dedup_s1_3 <- merge(test1_dedup, test3_dedup)
test_dedup_s1_4 <- merge(test1_dedup, test4_dedup)
test_dedup_s1_5 <- merge(test1_dedup, test5_dedup)
test_dedup_s1_6 <- merge(test1_dedup, test6_dedup)
test_dedup_s1_7 <- merge(test1_dedup, test7_dedup)
test_dedup_s1_8 <- merge(test1_dedup, test8_dedup)

test_dedup_s2_3 <- merge(test2_dedup, test3_dedup)
test_dedup_s2_4 <- merge(test2_dedup, test4_dedup)
test_dedup_s2_5 <- merge(test2_dedup, test5_dedup)
test_dedup_s2_6 <- merge(test2_dedup, test6_dedup)
test_dedup_s2_7 <- merge(test2_dedup, test7_dedup)
test_dedup_s2_8 <- merge(test2_dedup, test8_dedup)

test_dedup_s3_4 <- merge(test3_dedup, test4_dedup)
test_dedup_s3_5 <- merge(test3_dedup, test5_dedup)
test_dedup_s3_6 <- merge(test3_dedup, test6_dedup)
test_dedup_s3_7 <- merge(test3_dedup, test7_dedup)
test_dedup_s3_8 <- merge(test3_dedup, test8_dedup)

test_dedup_s4_5 <- merge(test4_dedup, test5_dedup)
test_dedup_s4_6 <- merge(test4_dedup, test6_dedup)
test_dedup_s4_7 <- merge(test4_dedup, test7_dedup)
test_dedup_s4_8 <- merge(test4_dedup, test8_dedup)

test_dedup_s5_6 <- merge(test5_dedup, test6_dedup)
test_dedup_s5_7 <- merge(test5_dedup, test7_dedup)
test_dedup_s5_8 <- merge(test5_dedup, test8_dedup)

test_dedup_s6_7 <- merge(test6_dedup, test7_dedup)
test_dedup_s6_8 <- merge(test6_dedup, test8_dedup)

test_dedup_s7_8 <- merge(test7_dedup, test8_dedup)


test_s1_dedup_all8 <- test1_dedup[which(!(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_2['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_3['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_4['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_5['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_6['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_7['X.Uploaded_variation'])) &
                                          !(unlist(test1_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_8['X.Uploaded_variation']))),]

test_s2_dedup_all8 <- test2_dedup[which(!(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_2['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_3['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_4['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_5['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_6['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_7['X.Uploaded_variation'])) &
                                          !(unlist(test2_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_8['X.Uploaded_variation']))),]

test_s3_dedup_all8 <- test3_dedup[which(!(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_3['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_3['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_4['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_5['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_6['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_7['X.Uploaded_variation'])) &
                                          !(unlist(test3_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_8['X.Uploaded_variation']))),]

test_s4_dedup_all8 <- test4_dedup[which(!(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_4['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_4['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_4['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_5['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_6['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_7['X.Uploaded_variation'])) &
                                          !(unlist(test4_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_8['X.Uploaded_variation']))),]

test_s5_dedup_all8 <- test5_dedup[which(!(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_5['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_5['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_5['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_5['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_6['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_7['X.Uploaded_variation'])) &
                                          !(unlist(test5_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_8['X.Uploaded_variation']))),]

test_s6_dedup_all8 <- test6_dedup[which(!(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_6['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_6['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_6['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_6['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_6['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s6_7['X.Uploaded_variation'])) &
                                          !(unlist(test6_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s6_8['X.Uploaded_variation']))),]

test_s7_dedup_all8 <- test7_dedup[which(!(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s6_7['X.Uploaded_variation'])) &
                                          !(unlist(test7_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s7_8['X.Uploaded_variation']))),]

test_s8_dedup_all8 <- test8_dedup[which(!(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_all_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s1_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s2_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s3_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s4_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s5_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s6_8['X.Uploaded_variation'])) &
                                          !(unlist(test8_dedup['X.Uploaded_variation']) %in% unlist(test_dedup_s7_8['X.Uploaded_variation']))),]

# save(test_s1_dedup_all8, test_s2_dedup_all8, test_s3_dedup_all8, test_s4_dedup_all8, test_s5_dedup_all8,
#      test_s6_dedup_all8, test_s7_dedup_all8, test_s8_dedup_all8, test_all_8, file = "mutations_all8_unique.RData")

# load(file = "mutations_all8_unique.RData")

## Load package needed to identify potentially impactful mutations
library(dndscv)
library(data.table)

getDatFormat <- function(data) {
  res <- str_match(data, "([^_]*)_(.*?)_(.*?)/(.*?)(/|$)")
}

## Get data in right format
muts_VU_all8 <- as.data.frame(do.call(rbind, lapply(test_s1_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_VU_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_VU_all8$SampleID <- "VU"

muts_MGH_all8 <- as.data.frame(do.call(rbind, lapply(test_s2_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_MGH_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_MGH_all8$SampleID <- "MGH"

muts_BR1_all8 <- as.data.frame(do.call(rbind, lapply(test_s3_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_BR1_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_BR1_all8$SampleID <- "BR1"

muts_DS3_all8 <- as.data.frame(do.call(rbind, lapply(test_s4_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS3_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS3_all8$SampleID <- "DS3"

muts_DS6_all8 <- as.data.frame(do.call(rbind, lapply(test_s5_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS6_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS6_all8$SampleID <- "DS6"

muts_DS7_all8 <- as.data.frame(do.call(rbind, lapply(test_s6_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS7_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS7_all8$SampleID <- "DS7"

muts_DS8_all8 <- as.data.frame(do.call(rbind, lapply(test_s7_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS8_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS8_all8$SampleID <- "DS8"

muts_DS9_all8 <- as.data.frame(do.call(rbind, lapply(test_s8_dedup_all8$X.Uploaded_variation, getDatFormat)))[,1:5]
names(muts_DS9_all8) <- c("SampleID", "chr", "pos", "ref", "mut") 
muts_DS9_all8$SampleID <- "DS9"

## Put in common dataframe
muts_all8 <- do.call("rbind", list(muts_VU_all8, muts_MGH_all8,
                                   muts_BR1_all8, muts_DS3_all8, 
                                   muts_DS6_all8, muts_DS7_all8, 
                                   muts_DS8_all8, muts_DS9_all8))

## Run the analysis
dndsout_all8 = dndscv(muts_all8, refdb="~/git/GES_2020/WES/RefCDS_human_GRCh38.p12.rda", 
                      cv=NULL, max_muts_per_gene_per_sample = 5)
all8_keymuts <- subset(dndsout_all8$sel_cv, pglobal_cv < 0.05)
impGenes_all8 <- subset(dndsout_all8$annotmuts, gene %in% all8_keymuts$gene_name)

## Change dNdScv categorization to more specific scheme (for reviewers)
impGenes_all8$impact[impGenes_all8$impact == "no-SNV" & !is.na(impGenes_all8$ntchange)] <- sapply(str_split(impGenes_all8$ntchange[impGenes_all8$impact == "no-SNV" & !is.na(impGenes_all8$ntchange)], "-",3), function(x) x[3])
impGenes_all8 <- subset(impGenes_all8, impact != "no-SNV")
## Remove gene with mutations in 2+ samples (MUC3A, GOLGA6L7, MADCAM1, GOLGA8F, CHD4, GOLGA6L4, IGFN1) - for reviewers
muts_remove <- c("MUC3A","GOLGA6L7","MADCAM1","GOLGA8F","CHD4","GOLGA6L4","IGFN1")
impGenes_all8 <- subset(impGenes_all8, !(gene %in% muts_remove))

## hm gene order; rev(unique(impGenes_all8$gene))

## Plot cell line version cohort (of the mutations identified from joint analysis)
CLV <- c("VU", "MGH", "BR1")
impGenes_CLV <- subset(impGenes_all8, sampleID %in% CLV)

## Setting the color setup for plots below
matches <- c("delfrshift" = rainbow(7)[1], "delinframe" = rainbow(7)[2], "Essential_Splice" = rainbow(7)[3],
             "insfrshift" = rainbow(7)[4], "Missense" = rainbow(7)[5], "Nonsense" = rainbow(7)[6],
             "Synonymous" = rainbow(7)[7], "None" = "white")
names <- c("delfrshift", "delinframe" , "Essential_Splice", "insfrshift", 
           "Missense", "Nonsense", "Synonymous", "None")
labels <- c("Frameshift Deletion", "Inframe Deletion", "Splice Site", "Frameshift Insertion", 
            "Missense", "Nonsense", "Synonymous", "None")

## Creating dataframe-mutation grid (for heatmap plotting)
fillInDF_CLV_dNdS <- expand.grid(unique(impGenes_all8$gene), CLV)
names(fillInDF_CLV_dNdS) <- c("gene", "sampleID")
fillInDF_CLV_dNdS$interaction <- paste0(fillInDF_CLV_dNdS$gene, ":", fillInDF_CLV_dNdS$sampleID)

## Creating compiled dataframe of mutations, amino acids, and impact
impGenes_CLV$interaction <- paste0(impGenes_CLV$gene, ":", impGenes_CLV$sampleID)

## Compile new dataframe for each mutated gene, CLV, and impact
Genes_CLV_dNdS <- merge(fillInDF_CLV_dNdS, impGenes_CLV, c("interaction"), all.x = TRUE)
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
CG_CLV_dNdS_geneCount <- as.data.frame(lapply(CG_CLV_dNdS_geneCount, unlist))


## Plot vertical barplot for each gene (gene and impact classified)
ggplot() +
  geom_bar(CG_CLV_dNdS_geneCount, stat = "identity", color = "black", 
           aes(x = gene, y = n, fill = as.factor(impact)), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  scale_y_continuous(breaks=seq(0, 5, by = 1)) + 
  scale_x_discrete(limits = rev(levels(CG_CLV_dNdS_factor_p$gene))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    # axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12),
    axis.title.x=element_blank(), axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), panel.border = element_blank()) +
  ggsave("FIG_3D_right.pdf", width = 12, height = 2)

plt_test <- ggplot() +
  geom_bar(CG_CLV_dNdS_geneCount, stat = "identity", color = "black", 
           aes(x = gene, y = n, fill = as.factor(impact)), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels)

plt_test_leg <- ggpubr::get_legend(plt_test)
ggpubr::as_ggplot(plt_test_leg) + ggsave("FIG_3D_legend.pdf")


## Plot horizontal barplot for each sample (sample and impact classified)
CLV_plot <- impGenes_CLV %>% count(sampleID, impact)
CLV_plot <- as.data.frame(lapply(CLV_plot, unlist))

ggplot() +
  geom_bar(CLV_plot, stat = "identity", color = "black", 
           aes(x = factor(sampleID, levels = c("BR1", "MGH", "VU")),
               y = n, fill = as.factor(impact)), size = 0.125, width = 0.99) +
  theme_bw() + labs(x = "Sample", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  scale_y_continuous(breaks=seq(0, 30, by = 10)) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12),
    panel.border = element_blank()) +
  ggsave("FIG_3D_top.pdf", width = 6, height = 3)

## Add character and numeric indecies for heatmap plotting
e <- rev(unique(c(as.character(CG_CLV_dNdS_factor_p$sampleID))))
f <- as.numeric(factor(CG_CLV_dNdS_factor_p$sampleID, levels=e))

## Plot heatmap of mutations and CLV (colored by impact and number of mutations)
### Separated bars within sample bars are if multiple mutation impacts in same CLV
CG_CLV_dNdS_factor_p <- as.data.frame(lapply(CG_CLV_dNdS_factor_p, unlist))

ggplot(CG_CLV_dNdS_factor_p, (aes(x = gene, y = f + shift,
                                  fill = as.factor(impact), alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = as.factor(impact), alpha = n)) + theme_bw() + 
  labs(x = "Genes", y = "Sample") +
  scale_alpha(range = c(1,1), limits = c(1,3), name = "Number of\nMutations") +
  scale_y_continuous(breaks = c(1,2,3), labels = e) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names,
                    labels = labels) +
  theme(
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=16),
    legend.title = element_text(size=12), axis.title=element_text(size=16)) + 
  ggsave("FIG_3D_middle.pdf", width = 6, height = 12)

## Export gene list to csv - Supplementary Table S2

### Sublines

sublines <- c("DS3", "DS6", "DS7", "DS8", "DS9")
impGenes_sublines <- subset(impGenes_all8, sampleID %in% sublines)

## Creating dataframe-mutation grid (for heatmap plotting)
fillInDF_sublines_dNdS <- expand.grid(unique(impGenes_all8$gene), sublines)
names(fillInDF_sublines_dNdS) <- c("gene", "sampleID")
fillInDF_sublines_dNdS$interaction <- paste0(fillInDF_sublines_dNdS$gene, ":", fillInDF_sublines_dNdS$sampleID)

## Creating compiled dataframe of mutations, amino acids, and impact
impGenes_sublines$interaction <- paste0(impGenes_sublines$gene, ":", impGenes_sublines$sampleID)

## Compile new dataframe for each mutated gene, sublines, and impact
Genes_sublines_dNdS <- merge(fillInDF_sublines_dNdS, impGenes_sublines, c("interaction"), all.x = TRUE)
Genes_sublines_dNdS <- Genes_sublines_dNdS[,c("gene.x", "sampleID.x", "impact")]
Genes_sublines_dNdS[is.na(Genes_sublines_dNdS)] <- "None"
names(Genes_sublines_dNdS) <- c("gene", "sampleID", "impact")

## Tally number of mutations for each sample-gene pair
CG_sublines_factor_dNdS <- Genes_sublines_dNdS %>% modify_if(is.character, as.factor)
CG_sublines_factor_dNdS_count <- CG_sublines_factor_dNdS %>% 
  group_by(gene, sampleID, impact) %>% 
  tally() 

## Add shift and height column (annotations to heatmap)
CG_sublines_dNdS_factor_p <- data.table(CG_sublines_factor_dNdS_count)
CG_sublines_dNdS_factor_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
CG_sublines_dNdS_factor_p[, height:=1/.N, by=list(gene, sampleID)]

## Add character and numeric indecies for heatmap plotting
k <- unique(c(as.character(CG_sublines_dNdS_factor_p$sampleID)))
l <- as.numeric(factor(CG_sublines_dNdS_factor_p$sampleID, levels=k)) 

## Create gene count dataframe (for vertical heatmap bar plot)
CG_sublines_dNdS_geneCount <- subset(CG_sublines_factor_dNdS, impact != "None") %>% 
  group_by(gene, impact) %>% 
  tally() 

### Change from factor to character (for plots)
CG_sublines_dNdS_geneCount$gene <- as.character(CG_sublines_dNdS_geneCount$gene) 
CG_sublines_dNdS_geneCount <- as.data.frame(lapply(CG_sublines_dNdS_geneCount, unlist))

## Plot vertical barplot for each gene (gene and impact classified)
ggplot() +
  geom_bar(CG_sublines_dNdS_geneCount, stat = "identity", color = "black", 
           aes(x = gene, y = n, fill = as.factor(impact)), size = 0.25) +
  theme_bw() + labs(x = "Gene", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  scale_y_continuous(breaks=seq(0, 5, by = 1)) + ylim(0,5) +
  scale_x_discrete(limits = rev(levels(CG_sublines_dNdS_factor_p$gene))) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12),
    axis.title.x=element_blank(), axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), panel.border = element_blank()) +
  ggsave("FIG_4D_right.pdf", width = 12, height = 2)

## Plot horizontal barplot for each sample (sample and impact classified)
sublines_plot <- impGenes_sublines %>% count(sampleID, impact)
sublines_plot <- as.data.frame(lapply(sublines_plot, unlist))

ggplot() +
  geom_bar(sublines_plot, stat = "identity", color = "black", 
           aes(x = factor(sampleID, levels = rev(sublines)),
               y = n, fill = as.factor(impact)), size = 0.125, width = 0.99) +
  theme_bw() + labs(x = "Sample", y = "Number of Mutations") +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  scale_y_continuous(breaks=seq(0, 30, by = 10)) + ylim(0,30) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_text(size=12), axis.title=element_text(size=12),
    panel.border = element_blank()) +
  ggsave("FIG_4D_top.pdf", width = 6, height = 3)

## Add character and numeric indecies for heatmap plotting
e <- rev(unique(c(as.character(CG_sublines_dNdS_factor_p$sampleID))))
f <- as.numeric(factor(CG_sublines_dNdS_factor_p$sampleID, levels=e))

## Plot heatmap of mutations and sublines (colored by impact and number of mutations)
### Separated bars within sample bars are if multiple mutation impacts in same sublines
CG_sublines_dNdS_factor_p <- as.data.frame(lapply(CG_sublines_dNdS_factor_p, unlist))

ggplot(CG_sublines_dNdS_factor_p, (aes(x = gene, y = f + shift,
                                       fill = as.factor(impact), alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = as.factor(impact), alpha = n)) + theme_bw() + 
  labs(x = "Genes", y = "Sample") +
  scale_alpha(range = c(1,1), limits = c(1,5), name = "Number of\nMutations") +
  scale_y_continuous(breaks = rev(seq(5)), labels = e) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names,
                    labels = labels) +
  theme(
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    legend.position = "none", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(size = 6, angle = 90, hjust = 0),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=16),
    legend.title = element_text(size=12), axis.title=element_text(size=16)) + 
  ggsave("FIG_4D_middle.pdf", width = 6, height = 12)


### Plot Cancer Associated Genes
## Cancer gene list
cancerGenes <- c("EGFR", "MET", "HER2", "BRAF", "NF1", "RAF1", "KRAS", "NRAS", "HRAS", "MAP2K2", 
                 "AKT1", "PIK3CA", "PIK3CB", "TSC1", "ALK", "APC", "EPHA2", "EPHA3", "EPHA4",
                 "ERBB4", "FGFR1", "ITK", "JAK2", "JAK3", "LRP1B", "LTK", "ROS1", "STK11", "TP53", "RB1", 
                 "ABCB5", "CFTR", "DACH1", "RELN", "CDKN2A", "DDR2", "DLEC1", "IRF1", "KEAP1", "MAP2K1",
                 "MAP3K8", "NRG1", "PPP2R1B", "PRKN", "PTEN", "RASSF1", "RET", "RIT1", "SLC22A18", "SMARCA4")

## Create empty dataframe
all8 <- c("BR1", "MGH", "VU", "DS3", "DS6", "DS7", "DS8", "DS9")
fillInDF_all8_CG <- expand.grid(cancerGenes, all8)
names(fillInDF_all8_CG) <- c("gene", "sampleID")
fillInDF_all8_CG$interaction <- paste0(fillInDF_all8_CG$gene, ":", fillInDF_all8_CG$sampleID)

## Subset data for cancer genes
all8_all_CG <- subset(dndsout_all8$annotmuts, gene %in% cancerGenes)
all8_all_CG$impact[all8_all_CG$impact == "no-SNV" & !is.na(all8_all_CG$ntchange)] <- sapply(str_split(all8_all_CG$ntchange[all8_all_CG$impact == "no-SNV" & !is.na(all8_all_CG$ntchange)], "-",3), function(x) x[3])
all8_all_CG <- subset(all8_all_CG, impact != "no-SNV")
all8_all_CG$interaction <- paste0(all8_all_CG$gene, ":", all8_all_CG$sampleID)

## Merge two dataframes together
cancerGenes_all8 <- merge(fillInDF_all8_CG, all8_all_CG, c("interaction"), all.x = TRUE)
cancerGenes_all8 <- cancerGenes_all8[,c("gene.x", "sampleID.x", "impact")]
cancerGenes_all8[is.na(cancerGenes_all8)] <- "None"
names(cancerGenes_all8) <- c("gene", "sampleID", "impact")

## Tally the total number of mutations for each interaction group
test_f <- cancerGenes_all8 %>% modify_if(is.character, as.factor)
test_f <- test_f %>% 
  group_by(gene, sampleID, impact) %>% 
  tally()

# Plot together
test_f_all_p <- data.table(test_f)
test_f_all_p[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
test_f_all_p[, height:=1/.N, by=list(gene, sampleID)]

## Create indecies
i <- unique(c(as.character(test_f_all_p$sampleID)))
j <- as.numeric(factor(test_f_all_p$sampleID, levels=i))

## Plot
ggplot(test_f_all_p, (aes(x = gene, y = j + shift,
                          fill = impact, alpha = n, height = height))) + 
  geom_tile(color = "black", aes(fill = impact, alpha = n)) + theme_bw() + 
  labs(x = "Unique Cancer Genes", y = "Sample") +
  scale_alpha(range = c(1,1)) +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8), 
                     labels = i) +
  coord_flip() +
  scale_fill_manual(values = matches, breaks = names, labels = labels) +
  theme(
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    legend.position = "right", axis.ticks = element_blank(),
    panel.border=element_rect(fill = NA, colour="white",size=0.2),
    axis.text.y = element_text(size = 12), legend.text = element_text(size = 16),
    plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=16),
    legend.title = element_text(size=12), axis.title=element_text(size=16)) + 
  ggsave("FIG_S4A.pdf", width = 12, height = 12)

# # Save data for GO semantic similarity analysis
# save(test_s1_dedup, test_s2_dedup, test_s3_dedup,
#      test_s4_dedup_noDS8, test_s5_dedup_noDS8, test_s6_dedup_noDS8, test_s8_dedup_noDS8,
#      test_s4_dedup, test_s5_dedup, test_s6_dedup, test_s7_dedup, test_s8_dedup,
#      file = "../GO/variants_byCohort.RData")