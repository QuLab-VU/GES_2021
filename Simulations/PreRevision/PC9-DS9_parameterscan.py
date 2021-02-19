# Load necessary packages
from pysb import *
import numpy as np
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import pandas as pd
import math

# Read in experimental cFP data for associated sublime
cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")
cFP_DS9_Erl = cFP_rates[(cFP_rates['Cell_Line'] == 'PC9-DS9')]

# Function to run a monoclonal growth model (MGM) simulation, 
# collect simulated DIP rates, and compare to experimental sublime
# cFP data
def dist_compare(dat, num_cells, div, dth):
    ## Initiate model with associated data and parameters
    Model()
    Monomer('Cell')
    Parameter('Cell_init', num_cells)
    Initial(Cell, Cell_init)
    Observable('Obs_Cell', Cell())
    Parameter('k_div', div)
    Parameter('k_dth', dth)
    Rule('Division', Cell() >> Cell() + Cell(), k_div)
    Rule('Death', Cell() >> None, k_dth)
    ## Run as many simulations as the number of colonies in associated subline cFP
    sim_dist = []
    for rep in range(len(dat)):
	 ## Setup time over which simulated colonies are allowed to grow
        t1 = np.linspace(0, 200, 201) # in hours
	 ## Simulate model over timespan - stochastic, so each result is slightly different
        sim1 = BngSimulator(model, tspan=t1, verbose=False)
	 ## Run simulation with ‘untreated’ division and death rate parameters
        x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                             'k_dth': 0.005 * np.log(2)},
                     verbose=False)
	 ## Collect the last time point measurement from ‘untreated’ colony growth 
	 ## phase of the simulation
        cell_tot = x1.observables["Obs_Cell"][-1]
	 ## Setup time over which to calculate colony DIP rate
        t2 = np.linspace(0, 225, 226)  # 8 days in drug
	 ## Run ‘treated’ phase of simulation with parameters defined in loop below
	 ## and with an initial conditions for each colony equal to the last time point
	 ## from the untreated phase
        x2 = sim1.run(tspan=t2, param_values={'k_div': model.parameters['k_div'].value,
                                             "k_dth": model.parameters['k_dth'].value,
                                             "Cell_init": cell_tot},
                     n_sim=1, verbose=False)
	 ## Collect DIP rates (slope of proliferation rate trajectory) and add to 
	 ## simulated DIP rate distribution count
        if cell_tot != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
    ## Compare the experimental and simulated distributions using the K-S test
    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    ## Return the p-value from the K-S test
    return p_val

# Enumerate the division and death rate range over which parameter scan 
# will be performed
divs_DS9 = np.linspace(0.018, 0.06, 42) * np.log(2)
dips_DS9 = np.linspace(0.0005, 0.0021, 16) * np.log(2)
# Initiate lists of each important parameter to be collected during the
# simulation
p_vals_DS9 = []
div_rates_DS9 = []
dth_rates_DS9 = []
dip_rates_DS9 = []
# Loop over each division and death (defined by division-death) rate,
# run the model, append data to the lists above
for dip in dips_DS9:
    for di in divs_DS9:
        dt = di-dip
        # print(dip, di, dt)
        p = dist_compare(cFP_DS9_Erl, 1, di, dt)
        p_vals_DS9.append(p)
        dip_rates_DS9.append(dip)
        div_rates_DS9.append(di)
        dth_rates_DS9.append(dt)
        print("run done")

# Create a dictionary with data collected above
dict = {'DIP rate': dip_rates_DS9, 'division rate': div_rates_DS9,
        'death rate': dth_rates_DS9, 'p-value': p_vals_DS9}

# Create a pandas data frame from dictionary and save as pickled object
df = pd.DataFrame(data=dict)
df.to_pickle('PC9-DS9_param-scan_tighterRange.pkl')