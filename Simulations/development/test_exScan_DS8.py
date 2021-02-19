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

def dist_compare_2S(samp, div, dth, len_ts):
    ## Pull out only rates from DS8
    dat = cFP_rates[cFP_rates['Cell_Line'] == samp]
    num_cells = 1
    Model()
    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(div))]})
    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(div))]
    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, div[i])) for i in range(len(div))]
    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, dth[i])) for i in range(len(div))]
    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(div))]
    ## Collect simulated DIP rate distributions in list
    sim_dist = []
    ## Run first part of simulation (pre-drug: set division and death rates)
    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim1.run(tspan=t1, param_values={'kdiv_0': 0.04 * np.log(2),
                                          'kdiv_1': 0.04 * np.log(2),
                                          'kdth_0': 0.005 * np.log(2),
                                          'kdth_1': 0.005 * np.log(2)},
                 n_runs=len(dat), verbose=False)
    ## Capture response trajectories from first part of simulation
    trajs = np.array(np.array([tr[:]["Cell_total"] for tr in np.array(x1.observables)]).T)
    ## Initiate second part of simulation (drug-treated) with last time point from first
    ## part of simulation
    cell_tot = trajs[-1]
    ## Split cells into two states (1 has 75% of cells, 2 has 25% cells)
    cell_tot_num = len(cell_tot)
    cell_tot_s = np.random.multinomial(cell_tot_num, [0.75,0.25])
    ### Two helpful objects - not used later in code though
    prob_list = ([0]*int(cell_tot_s[0])) + ([1]*int(cell_tot_s[1]))
    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    ### Randomly distribute cells according to proportional split between states above
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]
    ## Run second part of simulation and capture trajectories - state 1
    trajects = []
    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
        print("%d" % ind)
        cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, len_ts, len_ts+1)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist.append(slope)
        if cell > 50:
            trajects.append(np.log2(x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
    ## Run second part of simulation and capture trajectories - state 2
    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        print("%d" % ind)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, len_ts, len_ts+1)  # in drug
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist.append(slope)
        if cell > 50:
            trajects.append(np.log2(x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
    ## Capture p-value from K-S test comparison of experimental and simulated data
    bootresult = bootstrap(np.array(dat['DIP_Rate']), 100)

    ks_p_list = []
    ad_p_list = []
    for i in range(len(bootresult)):
        D,p = sp.ks_2samp(bootresult[i], sim_dist)
        st,c,s = sp.anderson_ksamp([bootresult[i], sim_dist])
        ks_p_list.append(p)
        ad_p_list.append(s)
    return trajects, sim_dist, ks_p_list, ad_p_list

# Run function for DS8 - PGM at specified example rates and capture key data
data_DS8 = dist_compare_2S(samp = "PC9-DS8", div = [0.022, 0.024],
                          dth = [0.0211, 0.0169], len_ts = 148)

## Output model trajectories
trajectories_DS8 = pd.DataFrame(data_DS8[0])
trajectories_DS8 = trajectories_DS8.transpose()
trajectories_DS8.to_csv('trajectories_DS8_G50_REVISION.csv')

## Output model distributions
distributions_DS8 = pd.DataFrame({'DS8': data_DS8[1]})
distributions_DS8.to_csv('distributions_DS8_G50_REVISION.csv')

## Output model p-values from KS test
KSbootstrap_DS8 = pd.DataFrame({'DS8': data_DS8[2]})
KSbootstrap_DS8.to_csv('KSbootstrap_DS8_G50_REVISION.csv')

## Output model p-values from AD test
ADbootstrap_DS8 = pd.DataFrame({'DS8': data_DS8[3]})
ADbootstrap_DS8.to_csv('ADbootstrap_DS8_G50_REVISION.csv')


