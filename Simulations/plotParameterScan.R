### Plot MGM and PGM parameter scans ###
library(ggplot2)

setwd('~/git/GES_2020/Simulations/')

# Load all of parameter scan data
tile <- read.csv('all_cellLine_tile.csv')

# Remove DS8 - only fit Polyclonal Growth Model (PGM - see below)
## Leaves only Monoclonal Growth Model (MGM) comparisons
tile_new <- subset(tile, !cell.line == "PC9.DS8")

# Set color and labels
cols_PS <- c("not.assigned" = "grey90", "PC9.DS1" = "coral", "PC9.DS3" = "brown",
             "PC9.DS4" = "deepskyblue", "PC9.DS6" = "deeppink", "PC9.DS7" = "darkorchid",
             "PC9.DS9" = "gold")
pops_PS <- c("not.assigned", "PC9.DS1", "PC9.DS3", "PC9.DS4", 
             "PC9.DS6","PC9.DS7", "PC9.DS9")
labs_PS <- c("Not Assigned", "DS1", "DS3", "DS4", "DS6", "DS7", "DS9")

# Plot all parameter scans (excluding DS8) - not in paper
ggplot(tile_new) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS1'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "coral") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS3'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "brown") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS4'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deepskyblue") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS6'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deeppink") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS7'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "darkorchid") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS9'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "gold") +
  theme_bw() + 
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value") +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("parameterScan_all.pdf", width = 10, height = 5)

# Plot subset of parameter scans (DS3 and DS4 - FIG. 5)
tile_forMainFig <- subset(tile, cell.line == c("PC9.DS3", "PC9.DS4")) 
cols <- c("PC9.DS3" = "brown", "PC9.DS4" = "deepskyblue")
pops <- c("PC9.DS3", "PC9.DS4")
labs <- c("DS3", "DS4")
ggplot(tile_forMainFig) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS3'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value),
            color = "brown") + #, fill = "brown") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS4'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), 
            color = "deepskyblue") + #, fill = "deepskyblue") +
  theme_bw() + 
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlab("Divison Rate - Death Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_5D.pdf", width = 10, height = 5)

# Plot different subset of parameter scans (DS1/6/7/9 - FIG. S9)
tile_forSuppFig <- subset(tile, cell.line %in% c("PC9.DS1", "PC9.DS6", "PC9.DS7", "PC9.DS9")) 
cols <- c("PC9.DS1" = "coral", "PC9.DS6" = "deeppink", 
          "PC9.DS7" = "darkorchid", "PC9.DS9" = "gold")
pops <- c("PC9.DS1", "PC9.DS6", "PC9.DS7", "PC9.DS9")
labs <- c("DS1", "DS6", "DS7", "DS9")
ggplot(tile_forSuppFig) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS1'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "coral") + 
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS6'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deeppink") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS7'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "darkorchid") + 
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS9'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "gold") + 
  theme_bw() +
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlim(-0.0002, 0.0022) +
  xlab("Division Rate - Death Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_S9D.pdf", width = 7.5, height = 4.5)


# Plot Polyclonal Growth Model (PGM) parameter scan on DS8

## Read in DS8 PGM data
tile_DS8_wholeRange <- read.csv('DS8_twoState_tile_wholeRange.csv')

# Split data into two states
DS8_WR_state1 <- subset(tile_DS8_wholeRange, Cell.Line == 'PC9-DS8.1')
DS8_WR_state2 <- subset(tile_DS8_wholeRange, Cell.Line == 'PC9-DS8.2')

## Add parameter pair metadata (due to scan in 4D - two division and two
## death rates - this is used to handle parameter value overlaps)
DS8_WR_state1$ParamPair <- as.character(seq(1:nrow(DS8_WR_state1)))
DS8_WR_state2$ParamPair <- as.character(seq(1:nrow(DS8_WR_state2)))

## Bind dataframes together
DS8_WR <- rbind(DS8_WR_state1, DS8_WR_state2)

## Create modified dataframe that jitters data points
### Prevents parameter overlaps that lead to plotting issues
DS8_WR_df <- data.frame(
  p.value = DS8_WR$p.value,
  Cell.Line = DS8_WR$Cell.Line,
  DIP.Rate = jitter(DS8_WR$DIP.Rate,2),
  Division.Rate = jitter(DS8_WR$Division.Rate,2),
  Death.Rate = jitter(DS8_WR$Death.Rate,2),
  ParamPair = as.numeric(DS8_WR$ParamPair)
)

# New color-label pattern
cols_DS8 <- c("PC9-DS8.1" = "seagreen", "PC9-DS8.2" = "seagreen")
pops_DS8 <- c("PC9-DS8.1", "PC9-DS8.2")
labs_DS8 <- c("DS8 State 1", "DS8 State 2")

plt_DS8 <- ggplot(DS8_WR_df, aes(DIP.Rate, Division.Rate, color = Cell.Line, 
                                 alpha = p.value, fill = Cell.Line)) +#, fill = ParamPair)) +
  geom_point(aes(alpha = p.value), color = "seagreen", fill = "seagreen") +
  theme_bw()  + ylim(0.01, 0.045) + xlim(0, 0.008) + 
  scale_color_manual(values = cols_DS8,
                     breaks = pops_DS8,
                     labels = labs_DS8,
                     name = "Subline") +
  scale_fill_manual(values = cols_DS8,
                    breaks = pops_DS8,
                    labels = labs_DS8,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlab("Division Rate - Death Rate") + ylab("Division Rate") +  
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

plt_DS8 + ggsave("FIG_6G.pdf", width = 5, height = 3)