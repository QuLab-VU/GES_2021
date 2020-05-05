# Script to convert pickled simulation data to information
# plotted in the parameter scan - DS8 two-state model

# Load necessary packages
import pandas as pd
import numpy as np

# Read in pickled simulated datasets
df_DS1 = pd.read_pickle('PC9-DS1_param-scan_tighterRange.pkl')
df_DS3 = pd.read_pickle('PC9-DS3_param-scan_tighterRange.pkl')
df_DS4 = pd.read_pickle('PC9-DS4_param-scan_tighterRange.pkl')
df_DS6 = pd.read_pickle('PC9-DS6_param-scan_tighterRange.pkl')
df_DS7 = pd.read_pickle('PC9-DS7_param-scan_tighterRange.pkl')
df_DS9 = pd.read_pickle('PC9-DS9_param-scan_tighterRange.pkl')

# Clean simulated data frames for each sublime

df_DS1 = df_DS1.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS1['cell.line'] = 'PC9.DS1'

df_DS3 = df_DS3.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS3['cell.line'] = 'PC9.DS3'

df_DS4 = df_DS4.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS4['cell.line'] = 'PC9.DS4'

df_DS6 = df_DS6.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS6['cell.line'] = 'PC9.DS6'

df_DS7 = df_DS7.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS7['cell.line'] = 'PC9.DS7'

df_DS9 = df_DS9.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS9['cell.line'] = 'PC9.DS9'

# Compile into common data frame for all sublines
all_dfs = [df_DS1, df_DS3, df_DS4, df_DS6, df_DS7, df_DS8, df_DS9]
df_all = pd.concat(all_dfs)

# Function to label cell line if achieve a certain p-value threshold
def label_cell_line(df):
    if df['p-value'] > 0.1:
        return df['cell.line']
    else:
        return "not.assigned"

# Add new column (corresponding to modified label - see above function)
df_all['cell.line.new'] = df_all.apply(lambda df_all: label_cell_line(df_all), axis = 1)

# Rename columns and save dataset
df_all.columns = ['DIP', 'death', 'division', 'pval', 'cellLine', 'cellLineNew']
df_all.to_csv('all_cellLine_tile.csv')