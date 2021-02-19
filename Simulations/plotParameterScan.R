library(ggplot2)
library(ggpubr)

setwd('~/git/GES_2020/Simulations/')

# Load all of parameter scan data
tile_DS1 <- read.csv('DS1_expansionTest_tile_lowVal.csv')
tile_DS3 <- read.csv('DS3_expansionTest_tile_lowVal.csv')
tile_DS4 <- read.csv('DS4_expansionTest_tile_lowVal.csv')
tile_DS6 <- read.csv('DS6_expansionTest_tile_lowVal.csv')
tile_DS7 <- read.csv('DS7_expansionTest_tile_lowVal.csv')
tile_DS9 <- read.csv('DS9_expansionTest_tile_lowVal.csv')

tile_all <- rbind(tile_DS1, tile_DS3, tile_DS4,
                  tile_DS6, tile_DS7, tile_DS9)

# Set color and labels
cols_PS <- c("not.assigned" = "white", "PC9.DS1" = "coral", "PC9.DS3" = "brown",
             "PC9.DS4" = "deepskyblue", "PC9.DS6" = "deeppink", "PC9.DS7" = "darkorchid",
             "PC9.DS9" = "gold")
pops_PS <- c("not.assigned", "PC9.DS1", "PC9.DS3", "PC9.DS4", 
             "PC9.DS6", "PC9.DS7", "PC9.DS9")
labs_PS <- c("Not Assigned", "DS1", "DS3", "DS4", "DS6", "DS7", "DS9")

# ggplot(tile_all) +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS1'),
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS3'),
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS4'),
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS6'),
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS7'),
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS9'), 
#             aes(DIP.rate, division.rate, fill = cell.line.new),
#             alpha = 0.7, color = "black") +
#   theme_bw() + ylim(0.01, 0.05) + xlim(-0.0020, 0.0040) +
#   scale_fill_manual(values = cols_PS,
#                     breaks = pops_PS,
#                     labels = labs_PS) +
#   xlab("DIP Rate") + ylab("Division Rate") +
#   theme(legend.text = element_text(size = 14), legend.position = "none", 
#         legend.title = element_text(size=14), 
#         axis.title=element_text(size=14), axis.text = element_text(size = 14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggsave("ParamScan_ADtestBootstrapped.pdf", width = 6, height = 4)


plt_scan_DS3_4 <- ggplot(tile_all) +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS3'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS4'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  theme_bw() + ylim(0.01, 0.05) + xlim(-0.0020, 0.0040) +
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS) +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 14), legend.position = "none", 
        legend.title = element_text(size=14), 
        axis.title=element_text(size=14), axis.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plt_scan_DS3_4 + ggsave("FIG_5D.pdf", width = 6, height = 4)
# leg_scan_DS3_4 <- get_legend(plt_scan_DS3_4)
# as_ggplot(leg_scan_DS3_4) + ggsave("legend_FIG_5D.pdf", width = 3, height = 1)

plt_scan_DS1_6_7_9 <- ggplot(tile_all) +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS1'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS6'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS7'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS9'), 
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  theme_bw() + ylim(0.01, 0.05) + xlim(-0.002, 0.0040) +
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS) +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 14), legend.position = "none", 
        legend.title = element_text(size=14), 
        axis.title=element_text(size=14), axis.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plt_scan_DS1_6_7_9 + ggsave("FIG_S12D.pdf", width = 6, height = 4)
# leg_scan_DS1_6_7_9 <- get_legend(plt_scan_DS1_6_7_9)
# as_ggplot(leg_scan_DS1_6_7_9) + ggsave("legend_FIG_S10D.pdf", width = 3, height = 1)


## Read in DS8 PGM data
tile_DS8_expandedRange <- read.csv('DS8_twoState_tile_expandedRange_lowVal_forPaper.csv')

# Split data into two states
DS8_WR_state1 <- subset(tile_DS8_expandedRange, Cell.Line == 'PC9-DS8.1')
DS8_WR_state2 <- subset(tile_DS8_expandedRange, Cell.Line == 'PC9-DS8.2')

## Add parameter pair metadata (due to scan in 4D - two division and two
## death rates - this is used to handle parameter value overlaps)
DS8_WR_state1$Param.Pair <- as.character(seq(1:nrow(DS8_WR_state1)))
DS8_WR_state2$Param.Pair <- as.character(seq(1:nrow(DS8_WR_state2)))

## Bind dataframes together
DS8_WR <- rbind(DS8_WR_state1, DS8_WR_state2)

## Create modified dataframe that jitters data points
### Prevents parameter overlaps that lead to plotting issues
DS8_WR_df <- data.frame(
  Cell.Line = DS8_WR$Cell.Line,
  Param.Pair = as.numeric(DS8_WR$Param.Pair),
  DIP.Rate = jitter(DS8_WR$DIP.Rate,2),
  Division.Rate = jitter(DS8_WR$Division.Rate,2),
  Death.Rate = jitter(DS8_WR$Death.Rate,2),
  KS.value = DS8_WR$KS.val,
  AD.value = DS8_WR$AD.val
)

DS8_WR_df <- subset(DS8_WR_df, AD.value > 0.05)

# New color-label pattern
cols_DS8 <- c("PC9-DS8.1" = "seagreen", "PC9-DS8.2" = "seagreen")
pops_DS8 <- c("PC9-DS8.1", "PC9-DS8.2")
labs_DS8 <- c("DS8 State 1", "DS8 State 2")

plt_DS8 <- ggplot(DS8_WR_df, aes(DIP.Rate, Division.Rate)) + #, fill = Cell.Line, group = as.factor(Param.Pair))) +
  geom_point(aes(color = Cell.Line)) +
  # geom_point(shape = 21, aes(color = as.factor(Param.Pair))) +
  # geom_line() +
  theme_bw() +
  scale_color_manual(values = cols_DS8,
                     breaks = pops_DS8,
                     labels = labs_DS8,
                     name = "Subline") +
  xlim(0, 0.008) + ylim(0, 0.08) +
  xlab("Division Rate - Death Rate") + ylab("Division Rate") +  
  theme(legend.text = element_text(size = 14), legend.position = "none",
        axis.text=element_text(size=14),
        legend.title = element_text(size=14), axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plt_DS8 + ggsave("FIG_5H.svg", 
                 width = 6, height = 4)
