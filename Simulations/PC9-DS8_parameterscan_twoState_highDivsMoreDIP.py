# Load necessary packages
from pysb import *
import numpy as np
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import pandas as pd
import math

# Load cFP experimental data
cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")

# Function to run two-state model and compare to DS8 cFP
def dist_compare(div, dth):
    ## Keep only DS8 data (for later comparison)
    dat = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS8']
    ## Define parameter values
    num_cells = 1
    kdiv = div
    kdth = dth
    ## Initiate model with parameters
    Model()
    ## Cell monomer is divided into two types - “states” for later
    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})
    ## Initial conditions defined for each “state”
    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]
    ## Each state has corresponding division and death equations
    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]
    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]
    ## Both a total cell count observable and each state defined
    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]
    
    ## Create empty list to add DIP rates into
    sim_dist = []
    ## Define time over which colonies proliferate in ‘untreated’ phase
    t1 = np.linspace(0, 200, 201) # hours
    ## Simulate ‘untreated’ phase of model
    sim1 = BngSimulator(model, tspan=t1, verbose=False)
    ## Run simulation with ‘untreated’ division and death rate parameters
    x1 = sim1.run(tspan=t1, param_values={'kdiv_0': 0.04 * np.log(2),
					  'kdiv_1': 0.04 * np.log(2),
                                          'kdth_0': 0.005 * np.log(2),
					  'kdth_1': 0.005 * np.log(2)},
                 n_runs=len(dat), verbose=False)
    ## Collect untreated trajectories 
    trajs = np.array(np.array([tr[:]["Cell_total"] for tr in np.array(x1.observables)]).T)
    ## Collect last time point from untreated trajectories
    cell_tot = trajs[-1]
    ## Multinomially distribute ‘colonies’ into two states with corresponding probabilities
    cell_tot_num = len(cell_tot)
    cell_tot_s = np.random.multinomial(cell_tot_num, [3/4.,1/4.])
    ## For each cell, provide the initial condition of the ‘treated’ phase
    ## of the simulation based on the state it corresponds to
    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]
    ## Loop over each colony in state 1 and run ’treated’ phase of simulation
    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
	 ### Keeping count of simulations
        print("%d" % ind)
	 ### Collect colony initial condition 
        cell_pop_1n = cell_pop1[ind]
	 ### Define timespan for ‘treated’ phase of simulation
        t2 = np.linspace(0, 225, 226)
	 ### Run simulation with associated parameters (see above)
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
	 ### Calculate DIP rates and collect into list
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
    ## Loop over each colony in state 2 and run ’treated’ phase of simulation
    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        print("%d" % ind)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, 225, 226)
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
    ## Compare simulated and experimental distributions using K-S test
    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    ## Return p-value of statistical test
    return p_val

# Creating range of division and dip rates to test over
di1 = np.linspace(0.036,0.040, 5)
di2 = np.linspace(0.036,0.040, 5)
dip1 = np.linspace(0.0005, 0.0015, 6)
dip2 = np.linspace(0.0055, 0.0075, 11)

# Create list of relevant parameters to be looped over for scan
ls_div = []
ls_dip = []
ls_dth = []
p_vals_DS8 = []
z = 0 # counter
# Run simulation for each division and death rate combination (two states - 4 dims)
for d1,d1_val in enumerate(di1):
    for d2, d2_val in enumerate(di2):
        for p1, p1_val in enumerate(dip1):
            for p2, p2_val in enumerate(dip2):
                ls_div.append([d1_val, d2_val])
                ls_dip.append([p1_val,p2_val])
                ls_dth.append([d1_val-p1_val, d2_val-p2_val])
                p = dist_compare(div = [d1_val, d2_val], dth = [d1_val-p1_val, d2_val-p2_val])
                p_vals_DS8.append(p)
                z = z + 1
                print(z)

# Create dictionary with data lists
dict = {'DIP rate': ls_dip, 'division rate': ls_div,
        'death rate': ls_dth, 'p-value': p_vals_DS8}

# Compile into pandas data frame and save as pickled object
df = pd.DataFrame(data=dict)
df.to_pickle('PC9-DS8_param-scan_twoState_highDivsMoreDIP.pkl')