from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import math
from astropy.stats import bootstrap

cFP_rates = pd.read_csv("/data/lola/hayforc/ParamScan/cFP_rates_VUlines.csv")

cFP_DS4_Erl = cFP_rates[(cFP_rates['Cell_Line'] == 'PC9-DS4')]

def dist_compare(dat, num_cells, div, dth):

    Model()
    Monomer('Cell')
    Parameter('Cell_init', num_cells)
    Initial(Cell, Cell_init)
    Observable('Obs_Cell', Cell())
    Parameter('k_div', div)
    Parameter('k_dth', dth)
    Rule('Division', Cell() >> Cell() + Cell(), k_div)
    Rule('Death', Cell() >> None, k_dth)
    sim_dist = []
    for rep in range(len(dat)):
        t1 = np.linspace(0, 200, 201) #hrs
        sim1 = BngSimulator(model, tspan=t1, verbose=False)
        x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                             'k_dth': 0.005 * np.log(2)},
                      verbose=False)
        cell_tot = x1.observables["Obs_Cell"][-1]
        t2 = np.linspace(0, 136, 137)  # length of experiment
        x2 = sim1.run(tspan=t2, param_values={'k_div': model.parameters['k_div'].value,
                                             "k_dth": model.parameters['k_dth'].value,
                                             "Cell_init": cell_tot},
                      n_sim=1, verbose=False)
        
        if cell_tot != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
    
    bootresult = bootstrap(np.array(dat['DIP_Rate']), 100)

    ks_p_list = []
    ad_p_list = []
    for i in range(len(bootresult)):
        D,p = sp.ks_2samp(bootresult[i], sim_dist)
        st,c,s = sp.anderson_ksamp([bootresult[i], sim_dist])
        ks_p_list.append(p)
        ad_p_list.append(s)

    return sim_dist, ks_p_list, ad_p_list

divs_DS4 = np.linspace(0.01, 0.08, 71)
dips_DS4 = np.linspace(0.0025, 0.0045, 21)
# divs_DS4 = np.linspace(0.034, 0.036, 3)
# dips_DS4 = np.linspace(0.0014, 0.0016, 3)

sim_dist_list = []
ks_ps = []
ad_ps = []
div_rates_DS4 = []
dth_rates_DS4 = []
dip_rates_DS4 = []

for dip in dips_DS4:
    for di in divs_DS4:
        dt = di-dip
        print(dip, di, dt)
        sim, ks, ad = dist_compare(cFP_DS4_Erl, 1, di, dt)
        sim_dist_list.append(sim)
        ks_ps.append(ks)
        ad_ps.append(ad)
        dip_rates_DS4.append(dip)
        div_rates_DS4.append(di)
        dth_rates_DS4.append(dt)
        print("run done")

dict = {'sim DIPs': sim_dist_list, 'DIP rate': dip_rates_DS4, 
        'division rate': div_rates_DS4, 'death rate': dth_rates_DS4, 
        'KS p-value': ks_ps, 'AD p-value': ad_ps}

df = pd.DataFrame(data=dict)

df.to_pickle('PC9-DS4_param-scan_Expansion.pkl')
