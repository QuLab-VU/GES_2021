#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import math
from astropy.stats import bootstrap

cFP_rates = pd.read_csv("/data/lola/hayforc/ParamScan/DSlines_expanded/bestRun/cFP_rates_VUlines.csv")



# In[3]:


def dist_compare(samp, div, dth):
    dat = cFP_rates[cFP_rates['Cell_Line'] == samp]
    num_cells = 1
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
    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)
    x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                          'k_dth': 0.005 * np.log(2)},
                  n_runs=len(dat), verbose=False)
    ## Capture response trajectories from first part of simulation
    trajs = np.array(np.array([tr[:]["Obs_Cell"] for tr in np.array(x1.observables)]).T)
    ## Initiate second part of simulation (drug-treated) with last time point from first
    ## part of simulation
    cell_pop = trajs[-1]
    ## Run second part of simulation and capture trajectories
    trajects = []
    for ind, cell in enumerate(cell_pop):
       	print("%d" % ind)
        print("%d" % cell)
       	t2 = np.linspace(0, 118, 119)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[0]: cell},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
               	x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
               	sim_dist.append(slope)
        if cell > 50:
            trajects.append(np.log2(x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
    ## Capture p-value from statistical test comparisons of experimental and simulated data
    bootresult = bootstrap(np.array(dat['DIP_Rate']), 100)

    ks_p_list = []
    ad_p_list = []
    for i in range(len(bootresult)):
        D,p = sp.ks_2samp(bootresult[i], sim_dist)
        st,c,s = sp.anderson_ksamp([bootresult[i], sim_dist])
        ks_p_list.append(p)
        ad_p_list.append(s)
    return trajects, sim_dist, ks_p_list, ad_p_list


# In[4]:


# Run function for sublines (except DS8 - PGM) at specified example rates
# and capture key data
DSs = ['PC9-DS1', 'PC9-DS3', 'PC9-DS4', 'PC9-DS6', 'PC9-DS7', 'PC9-DS9']
divs_a = [0.0200, 0.0150, 0.0250, 0.0200, 0.0170, 0.0150]
dths_a = [0.0190, 0.0155, 0.0214, 0.0196, 0.0150, 0.0138]
data_a = []
for i in range(len(DSs)):
    data_a.append(dist_compare(DSs[i], divs_a[i], dths_a[i]))


# In[8]:


## Output model trajectories
trajectories_DS1 = pd.DataFrame(data_a[0][0])
trajectories_DS1 = trajectories_DS1.transpose()
trajectories_DS1.to_csv('trajectories_DS1_G50_REVISION.csv')

trajectories_DS3 = pd.DataFrame(data_a[1][0])
trajectories_DS3 = trajectories_DS3.transpose()
trajectories_DS3.to_csv('trajectories_DS3_G50_REVISION.csv')

trajectories_DS4 = pd.DataFrame(data_a[2][0])
trajectories_DS4 = trajectories_DS4.transpose()
trajectories_DS4.to_csv('trajectories_DS4_G50_REVISION.csv')

trajectories_DS6 = pd.DataFrame(data_a[3][0])
trajectories_DS6 = trajectories_DS6.transpose()
trajectories_DS6.to_csv('trajectories_DS6_G50_REVISION.csv')

trajectories_DS7 = pd.DataFrame(data_a[4][0])
trajectories_DS7 = trajectories_DS7.transpose()
trajectories_DS7.to_csv('trajectories_DS7_G50_REVISION.csv')

trajectories_DS9 = pd.DataFrame(data_a[5][0])
trajectories_DS9 = trajectories_DS9.transpose()
trajectories_DS9.to_csv('trajectories_DS9_G50_REVISION.csv')



# In[9]:


## Output model distributions
distributions_0 = pd.DataFrame({'DS1': data_a[0][1]})
distributions_1 = pd.DataFrame({'DS3': data_a[1][1]})
distributions_2 = pd.DataFrame({'DS4': data_a[2][1]})
distributions_3 = pd.DataFrame({'DS6': data_a[3][1]})
distributions_4 = pd.DataFrame({'DS7': data_a[4][1]})
distributions_5 = pd.DataFrame({'DS9': data_a[5][1]})
distributions = pd.concat([distributions_0, distributions_1,
                           distributions_2, distributions_3,
                           distributions_4, distributions_5],
                          ignore_index=True, axis = 1)
distributions.columns = ['DS1', 'DS3', 'DS4', 'DS6', 'DS7', 'DS9']
distributions.to_csv('distributions_G50_REVISION.csv')



# In[22]:


## Output model p-values from KS test
KSbootstrap_0 = pd.DataFrame({'DS1': data_a[0][2]})
KSbootstrap_1 = pd.DataFrame({'DS3': data_a[1][2]})
KSbootstrap_2 = pd.DataFrame({'DS4': data_a[2][2]})
KSbootstrap_3 = pd.DataFrame({'DS6': data_a[3][2]})
KSbootstrap_4 = pd.DataFrame({'DS7': data_a[4][2]})
KSbootstrap_5 = pd.DataFrame({'DS9': data_a[5][2]})
KSbootstrap = pd.concat([KSbootstrap_0, KSbootstrap_1,
                         KSbootstrap_2, KSbootstrap_3,
                         KSbootstrap_4, KSbootstrap_5],
                        ignore_index=True, axis = 1)
KSbootstrap.columns = ['DS1', 'DS3', 'DS4', 'DS6', 'DS7', 'DS9']
KSbootstrap.to_csv('KSbootstrap_G50_REVISION.csv')


# In[23]:


## Output model p-values from KS test
ADbootstrap_0 = pd.DataFrame({'DS1': data_a[0][3]})
ADbootstrap_1 = pd.DataFrame({'DS3': data_a[1][3]})
ADbootstrap_2 = pd.DataFrame({'DS4': data_a[2][3]})
ADbootstrap_3 = pd.DataFrame({'DS6': data_a[3][3]})
ADbootstrap_4 = pd.DataFrame({'DS7': data_a[4][3]})
ADbootstrap_5 = pd.DataFrame({'DS9': data_a[5][3]})
ADbootstrap = pd.concat([ADbootstrap_0, ADbootstrap_1,
                         ADbootstrap_2, ADbootstrap_3,
                         ADbootstrap_4, ADbootstrap_5],
                        ignore_index=True, axis = 1)
ADbootstrap.columns = ['DS1', 'DS3', 'DS4', 'DS6', 'DS7', 'DS9']
ADbootstrap.to_csv('ADbootstrap_G50_REVISION.csv')


# In[24]:


ADbootstrap


# In[ ]:




