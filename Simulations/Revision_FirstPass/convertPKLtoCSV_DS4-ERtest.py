# Script to convert pickled simulation data to information
# plotted in the parameter scan - DS8 two-state model

# Load necessary packages
import pandas as pd
import numpy as np

# Read in pickled simulated datasets
df_D1LD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs1-LeftDIP.pkl')
df_D2LD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs2-LeftDIP.pkl')
df_D3LD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs3-LeftDIP.pkl')
df_D1RD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs1-RightDIP.pkl')
df_D2RD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs2-RightDIP.pkl')
df_D3RD = pd.read_pickle('PC9-DS4_param-scan_twoState-ER_Divs3-RightDIP.pkl')

# Clean simulated data frames for each sublime

df_D1LD = df_D1LD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D1LD['cell.line'] = 'PC9.D1LD'

df_D2LD = df_D2LD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D2LD['cell.line'] = 'PC9.D2LD'

df_D3LD = df_D3LD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D3LD['cell.line'] = 'PC9.D3LD'

df_D1RD = df_D1RD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D1RD['cell.line'] = 'PC9.D1RD'

df_D2RD = df_D2RD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D2RD['cell.line'] = 'PC9.D2RD'

df_D3RD = df_D3RD.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_D3RD['cell.line'] = 'PC9.D3RD'

# Compile into common data frame for all sublines
all_dfs = [df_D1LD, df_D2LD, df_D3LD, df_D1RD, df_D2RD, df_DS8, df_D3RD]
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
df_all.to_csv(‘DS4_ERtest_tile.csv')