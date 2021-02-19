# Load packages and read in data
import numpy as np
import pandas as pd
import scipy as sp
# from sklearn.neighbors import KernelDensity

# Input .pkl file for specified dateset
## DS1 used here
pkl_file = 'PC9-DS1_param-scan_Expansion.pkl'
population = 'PC9.DS1'
cFP_rates = pd.read_csv("/Users/Corey/git/GES_2020/Simulations/cFP_rates_VUlines.csv")
dat_DS1 = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS1']
data = dat_DS1

# Run function to identify lower bounds of bootstrapped p-value
def analyze_model(pkl_file, population, data):
    ## Read in dataset
    df = pd.read_pickle(pkl_file)
    df['cell line'] = population

    ## Create pandas objects from simulated data parameters
    df_dip = pd.DataFrame(df['DIP rate'].values.tolist(), columns=['DIP'])
    df_div = pd.DataFrame(df['division rate'].values.tolist(), columns=['div'])
    df_dth = pd.DataFrame(df['death rate'].values.tolist(), columns=['dth'])
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

# Run function
DS1_E = analyze_model(pkl_file=pkl_file, population = population, data = dat_DS1)


# Save dataset as CSV (for plotting)
DS1_E.columns = ['DIP rate', 'division rate', 'death rate', 'cell line', 'KS val', 'AD val']
DS1_E['cell.line.new'] = np.where(DS1_E['AD val'] > 0.05, population, 'not.assigned')
DS1_E = DS1_E[['cell line', 'division rate', 'death rate', 'DIP rate', 'KS val', 'AD val', 'cell.line.new']]
DS1_E.to_csv('DS1_expansionTest_tile_lowVal.csv')