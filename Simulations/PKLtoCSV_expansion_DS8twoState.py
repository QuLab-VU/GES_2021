# Script to convert pickled simulation data to information
# plotted in the parameter scan - DS8 two-state model

# Load necessary packages
import pandas as pd
import numpy as np

# Function to 
def analyze_model_DS8(pkl_file, population): #, data, num_params):
    ## Read in dataset
    df = pd.read_pickle(pkl_file)
    df['cell line'] = population

    ## Create pandas objects from simulated data parameters
    df_dip = pd.DataFrame(df['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
    df_div = pd.DataFrame(df['division rate'].values.tolist(), columns=['div1','div2'])
    df_dth = pd.DataFrame(df['death rate'].values.tolist(), columns=['dth1','dth2'])
    df_cellline = df['cell line']

    sim = np.array(df['sim DIPs'])

    KSval = []
    ADval = []


    for index, row in df.iterrows():
        KSval.append(np.mean(np.array(df['KS p-value'][index])) - 1*np.std(np.array(df['KS p-value'][index])))
        ADval.append(np.mean(np.array(df['AD p-value'][index])) - 1*np.std(np.array(df['AD p-value'][index])))

    df_c = pd.concat([df_dip.reset_index(drop=True),
                      df_div.reset_index(drop=True),
                      df_dth.reset_index(drop=True),
                      df_cellline.reset_index(drop=True)],
                     axis = 1)

    df_c['KS val'] = KSval
    df_c['AD val'] = ADval
    
    return(df_c)


# Load all subsets of model data
df_DS8_00_1 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs00_1_allDIP.pkl', population='PC9.DS8')
df_DS8_00_2 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs00_2_allDIP.pkl', population='PC9.DS8')
df_DS8_1 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs1_allDIP.pkl', population='PC9.DS8')
df_DS8_2 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs2_allDIP.pkl', population='PC9.DS8')
df_DS8_3 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs3_allDIP.pkl', population='PC9.DS8')
df_DS8_4 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs4_allDIP.pkl', population='PC9.DS8')
df_DS8_5 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs5_allDIP.pkl', population='PC9.DS8')
df_DS8_6 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs6_allDIP.pkl', population='PC9.DS8')
df_DS8_7 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs7_allDIP.pkl', population='PC9.DS8')
df_DS8_8 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs8_allDIP.pkl', population='PC9.DS8')
df_DS8_9 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs9_allDIP.pkl', population='PC9.DS8')
df_DS8_10 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs10_allDIP.pkl', population='PC9.DS8')
df_DS8_11 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs11_allDIP.pkl', population='PC9.DS8')
df_DS8_12 = analyze_model_DS8(pkl_file='PC9-DS8_param-scan_twoState_Divs12_allDIP.pkl', population='PC9.DS8')


# Compile into a single common data frame
DS8_all = pd.concat([df_DS8_00_1, df_DS8_00_2, df_DS8_1, df_DS8_2,
                    df_DS8_3, df_DS8_4, df_DS8_5, df_DS8_6,
                    df_DS8_7, df_DS8_8, df_DS8_9, df_DS8_10,
                    df_DS8_11, df_DS8_12], ignore_index = True)


# Annotate with modified identifier (ease of plotting)
DS8_all['cell line'] = np.where(DS8_all['AD val']>0.05, 'PC9-DS8', 'not.assigned')
DS8_all['param pair'] = range(DS8_all.shape[0])


# Create subsets of data frame (see below)
DS8_sig1 = DS8_all[['DIP1', 'div1', 'dth1', 'cell line', 'param pair', 'KS val', 'AD val']]
DS8_sig2 = DS8_all[['DIP2', 'div2', 'dth2', 'cell line', 'param pair', 'KS val', 'AD val']]


# Rename columns of data frame subsets
DS8_sig1.rename(columns={'DIP1': 'DIP Rate',
                         'div1': 'Division Rate',
                         'dth1': 'Death Rate'},
                 inplace=True)
DS8_sig2.rename(columns={'DIP2': 'DIP Rate',
                         'div2': 'Division Rate',
                         'dth2': 'Death Rate'},
                inplace=True)


# Separate subsets based on which DS8 state it resides in
DS8_sig1['Cell Line'] = np.where(DS8_sig1['cell line'] == "PC9-DS8", 'PC9-DS8.1', 'not.assigned')
DS8_sig2['Cell Line'] = np.where(DS8_sig2['cell line'] == "PC9-DS8", 'PC9-DS8.2', 'not.assigned')

# Concatenate datasets back together
DS8_sig_all = pd.concat([DS8_sig1, DS8_sig2])

# Save concatenated datasets into common CSV (for plotting)
DS8_sig_all.to_csv('DS8_twoState_tile_expandedRange_lowVal_forPaper.csv')