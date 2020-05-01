from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns
import math

cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")

cFP_DS6_Erl = cFP_rates[(cFP_rates['Cell_Line'] == 'PC9-DS6')]

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
        sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
        x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                             'k_dth': 0.005 * np.log(2)},
                     verbose=False)
        cell_tot = x1.observables["Obs_Cell"][-1]
        t2 = np.linspace(0, 225, 226)  # 8 days in drug
        x2 = sim1.run(tspan=t2, param_values={'k_div': model.parameters['k_div'].value,
                                             "k_dth": model.parameters['k_dth'].value,
                                             "Cell_init": cell_tot},
                     n_sim=1, verbose=False)

        if cell_tot != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
                # print "Normal", slope

    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)

    return p_val


divs_DS6 = np.linspace(0.018, 0.06, 42) * np.log(2)
dips_DS6 = np.linspace(-0.0007, 0.0013, 20) * np.log(2)
p_vals_DS6 = []
div_rates_DS6 = []
dth_rates_DS6 = []
dip_rates_DS6 = []
for dip in dips_DS6:
    for di in divs_DS6:
        dt = di-dip
        # print(dip, di, dt)
        p = dist_compare(cFP_DS6_Erl, 1, di, dt)
        p_vals_DS6.append(p)
        dip_rates_DS6.append(dip)
        div_rates_DS6.append(di)
        dth_rates_DS6.append(dt)
        print("run done")

dict = {'DIP rate': dip_rates_DS6, 'division rate': div_rates_DS6,
        'death rate': dth_rates_DS6, 'p-value': p_vals_DS6}

df = pd.DataFrame(data=dict)

df.to_pickle('PC9-DS6_param-scan_tighterRange.pkl')
