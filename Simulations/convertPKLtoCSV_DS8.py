# Script to convert pickled simulation data to information
# plotted in the parameter scan - DS8 two-state model

# Load necessary packages
import pandas as pd
import numpy as np

# Read in pickled simulated datasets
df_DS8_low = pd.read_pickle('PC9-DS8_param-scan_twoState_lowDivsMoreDIP.pkl')
df_DS8_lowMed = pd.read_pickle('PC9-DS8_param-scan_twoState_lowMedDivsMoreDIP.pkl')
df_DS8_medHigh = pd.read_pickle('PC9-DS8_param-scan_twoState_medHighDivsMoreDIP.pkl')
df_DS8_high = pd.read_pickle('PC9-DS8_param-scan_twoState_highDivsMoreDIP.pkl')

## Add cell line annotation to data
df_DS8_low['cell.line'] = 'PC9.DS8'
df_DS8_lowMed['cell.line'] = 'PC9.DS8'
df_DS8_medHigh['cell.line'] = 'PC9.DS8'
df_DS8_high['cell.line'] = 'PC9.DS8'

## Create pandas objects from simulated data parameters (separated by each
## parameter scan ranges (low, lowMed, medHigh, high)
### low
DS8_dipdf_low = pd.DataFrame(df_DS8_low['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_low = pd.DataFrame(df_DS8_low['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_low = pd.DataFrame(df_DS8_low['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_low = df_DS8_low[['p-value', 'cell.line']]
DS8_rest_low.reset_index(drop=True)

### lowMed
DS8_dipdf_lowMed = pd.DataFrame(df_DS8_lowMed['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_lowMed = pd.DataFrame(df_DS8_lowMed['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_lowMed = pd.DataFrame(df_DS8_lowMed['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_lowMed = df_DS8_lowMed[['p-value', 'cell.line']]
DS8_rest_lowMed.reset_index(drop=True)

### medHigh
DS8_dipdf_medHigh = pd.DataFrame(df_DS8_medHigh['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_medHigh = pd.DataFrame(df_DS8_medHigh['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_medHigh = pd.DataFrame(df_DS8_medHigh['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_medHigh = df_DS8_medHigh[['p-value', 'cell.line']]
DS8_rest_medHigh.reset_index(drop=True)

### high
DS8_dipdf_high = pd.DataFrame(df_DS8_high['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_high = pd.DataFrame(df_DS8_high['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_high = pd.DataFrame(df_DS8_high['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_high = df_DS8_high[['p-value', 'cell.line']]
DS8_rest_high.reset_index(drop=True)

## Concatenate each data object into common dataframe for each
## parameter scan region 
DS8_low = pd.concat([DS8_rest_low.reset_index(drop=True), 
                     DS8_dipdf_low.reset_index(drop=True), 
                     DS8_divdf_low.reset_index(drop=True),
                     DS8_dthdf_low], axis = 1)

DS8_lowMed = pd.concat([DS8_rest_lowMed.reset_index(drop=True), 
                     DS8_dipdf_lowMed.reset_index(drop=True), 
                     DS8_divdf_lowMed.reset_index(drop=True),
                     DS8_dthdf_lowMed], axis = 1)

DS8_medHigh = pd.concat([DS8_rest_medHigh.reset_index(drop=True), 
                     DS8_dipdf_medHigh.reset_index(drop=True), 
                     DS8_divdf_medHigh.reset_index(drop=True),
                     DS8_dthdf_medHigh], axis = 1)

DS8_high = pd.concat([DS8_rest_high.reset_index(drop=True), 
                     DS8_dipdf_high.reset_index(drop=True), 
                     DS8_divdf_high.reset_index(drop=True),
                     DS8_dthdf_high], axis = 1)

## Compile into a single common data frame
DS8_all = pd.concat([DS8_low, DS8_lowMed,
                      DS8_medHigh, DS8_high], ignore_index = True)
DS8 = DS8_all

## Annotate with modified identifier (ease of plotting)
DS8['cell.line'] = np.where(DS8['p-value']>0.1, 'PC9-DS8', 'not.assigned')
DS8['param.pair'] = range(DS8.shape[0])

## Create subsets of data frame (see below)
DS8_sig1 = DS8[['p-value', 'cell.line', 'DIP1', 'div1', 'dth1', 'param.pair']]
DS8_sig2 = DS8[['p-value', 'cell.line', 'DIP2', 'div2', 'dth2', 'param.pair']]

## Rename columns of data frame subsets
DS8_sig1.rename(columns={'DIP1': 'DIP Rate',
                         'div1': 'Division Rate',
                         'dth1': 'Death Rate'},
                 inplace=True)
DS8_sig2.rename(columns={'DIP2': 'DIP Rate',
                         'div2': 'Division Rate',
                         'dth2': 'Death Rate'},
                inplace=True)

## Separate subsets based on which DS8 state it resides in
DS8_sig1['Cell Line'] = np.where(DS8_sig1['cell.line'] == "PC9-DS8", 'PC9-DS8.1', 'not.assigned')
DS8_sig2['Cell Line'] = np.where(DS8_sig2['cell.line'] == "PC9-DS8", 'PC9-DS8.2', 'not.assigned')

## Concatenate datasets back together
DS8_sig_all = pd.concat([DS8_sig1, DS8_sig2])

## Save concatenated datasets into common CSV (for plotting)
DS8_sig_all.to_csv('DS8_twoState_tile_wholeRange.csv')