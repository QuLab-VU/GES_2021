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

def dist_compare(div, dth, prop1, prop2):
    dat = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS8']
    num_cells = 1
    kdiv = div
    kdth = dth

    Model()

    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})

    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]

    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]

    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]

    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]


    sim_dist = []

    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim1.run(tspan=t1, param_values={'kdiv_0': 0.04 * np.log(2),
					  'kdiv_1': 0.04 * np.log(2),
                                          'kdth_0': 0.005 * np.log(2),
					  'kdth_1': 0.005 * np.log(2)},
                 n_runs=len(dat), verbose=False)

    trajs = np.array(np.array([tr[:]["Cell_total"] for tr in np.array(x1.observables)]).T)

    cell_tot = trajs[-1]
    cell_tot_num = len(cell_tot)
    cell_tot_s = np.random.multinomial(cell_tot_num, [prop1,prop2])


    prob_list = ([0]*int(cell_tot_s[0])) + ([1]*int(cell_tot_s[1]))

    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]


    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
        # print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, 148, 149)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist.append(slope)

    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        # print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, 148, 149)  # in drug
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist.append(slope)

    bootresult = bootstrap(np.array(dat['DIP_Rate']), 100)

    ks_p_list = []
    ad_p_list = []
    for i in range(len(bootresult)):
        D,p = sp.ks_2samp(bootresult[i], sim_dist)
        st,c,s = sp.anderson_ksamp([bootresult[i], sim_dist])
        ks_p_list.append(p)
        ad_p_list.append(s)

    return sim_dist, ks_p_list, ad_p_list


di1 = np.linspace(0.066,0.070, 5)
di2 = np.linspace(0.066,0.070, 5)
dip1 = np.linspace(0.0005, 0.0015, 6)
dip2 = np.linspace(0.0055, 0.0075, 11)

sim_dist_list = []
ks_ps = []
ad_ps = []
ls_div = []
ls_dip = []
ls_dth = []

z = 0
for d1,d1_val in enumerate(di1):
    for d2, d2_val in enumerate(di2):
        for p1, p1_val in enumerate(dip1):
            for p2, p2_val in enumerate(dip2):
               	ls_div.append([d1_val, d2_val])
               	ls_dip.append([p1_val,p2_val])
               	ls_dth.append([d1_val-p1_val, d2_val-p2_val])
               	sim, ks, ad = dist_compare(div = [d1_val, d2_val], dth = [d1_val-p1_val, d2_val-p2_val], prop1 = 0.75, prop2 = 0.25)
               	sim_dist_list.append(sim)
                ks_ps.append(ks)
                ad_ps.append(ad)
               	z = z + 1
               	print(z)

dict = {'DIP rate': ls_dip, 'division rate': ls_div,
        'death rate': ls_dth, 'sim DIPs': sim_dist_list,
        'KS p-value': ks_ps, 'AD p-value': ad_ps}

df = pd.DataFrame(data=dict)

df.to_pickle('PC9-DS8_param-scan_twoState_Divs10_allDIP.pkl')
