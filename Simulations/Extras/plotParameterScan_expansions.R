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

ggplot(tile_all) +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS1'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS3'),
            aes(DIP.rate, division.rate, fill = cell.line.new),
            alpha = 0.7, color = "black") +
  geom_tile(data = subset(tile_all, cell.line.new == 'PC9.DS4'),
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
  theme_bw() + ylim(0.01, 0.05) + xlim(-0.0020, 0.0040) +
  scale_fill_manual(values = cols_PS,
                    breaks = pops_PS,
                    labels = labs_PS) +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 14), legend.position = "none", 
        legend.title = element_text(size=14), 
        axis.title=element_text(size=14), axis.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("ParamScan_ADtestBootstrapped.pdf", width = 6, height = 4)


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

plt_scan_DS3_4 + ggsave("ParamScan_DS3-DS4_ADtestBootstrapped.pdf", width = 6, height = 4)
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

plt_scan_DS1_6_7_9 + ggsave("ParamScan_DS1-DS6-DS7-DS9_ADtestBootstrapped.pdf", width = 6, height = 4)
# leg_scan_DS1_6_7_9 <- get_legend(plt_scan_DS1_6_7_9)
# as_ggplot(leg_scan_DS1_6_7_9) + ggsave("legend_FIG_S10D.pdf", width = 3, height = 1)
