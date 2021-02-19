library(ggplot2)
setwd('~/git/GES_2020/Simulations/')

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

plt_DS8 + ggsave("ParamScan_twoStateDS8_ADtestBootstrapped.svg", 
                 width = 6, height = 4)
