from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
#import seaborn as sns
import math
#import matplotlib



cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")


def dist_compare(div, dth):
    dat = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS8']
    num_cells = 1
    # kdiv = [0.026, 0.026]
    # kdth = [0.025, 0.020]
    kdiv = div
    kdth = dth



    Model()

    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})

    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]

    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]

    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]

    # [degrade(Cell(type=str(i)), Parameter('kdeg_%d' % i, kdeg[i])) for i in range(len(kdiv))]
    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]


    sim_dist = []
    # plt.figure(figsize=[4,3])

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
    cell_tot_s = np.random.multinomial(cell_tot_num, [3/4.,1/4.])


    prob_list = ([0]*int(cell_tot_s[0])) + ([1]*int(cell_tot_s[1]))

    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]


    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
        print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, 225, 226)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Cell_total"]), color = 'red', lw=1)

    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, 225, 226)  # in drug
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        # plt.plot(t1[-1] + x3.tout[0], np.log2(x3.observables["Cell_total"]), color = 'blue', lw=1)

    # plt.xlabel("Time (hours)")
    # plt.ylabel("Log2 Cell Count")

    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    mean_exp = np.mean(dat['DIP_Rate'])
    sd_exp = np.std(dat['DIP_Rate'])
    mean_sim = np.mean(sim_dist)
    sd_sim = np.std(sim_dist)

    return p_val
    # plt.figure(figsize=[4,3])
    # sns.distplot(sim_dist, kde = True, hist= False, color = "grey", kde_kws={"shade": True})
    # sns.distplot(dat['DIP_Rate'], kde = True, hist= False, color = "seagreen", kde_kws={"shade": True})
    # plt.xlabel("DIP Rate")
    # plt.ylabel("Density")
    # plt.xlim(-0.025, 0.025)
    # plt.ylim(0,200)
    # plt.legend(labels=['Simulated','Experimental'], loc = 'upper left')
    # plt.text(0.007, 200, "Mean (Exp)=%.3f" % round(mean_exp, 3), fontsize = 8)
    # plt.text(0.007, 180, "SD (Exp)=%.4f" % round(sd_exp, 4), fontsize = 8)
    # plt.text(0.007, 160, "Mean (Sim)=%.3f" % round(mean_sim, 3), fontsize = 8)
    # plt.text(0.007, 140, "SD (Sim)=%.4f" % round(sd_sim, 4), fontsize = 8)
    # plt.text(0.007, 150, "p=%.3f" % round(p_val, 3), fontsize = 12)
    # plt.title("DIP Rate Distribution")


di1 = np.linspace(0.026,0.030, 5)
di2 = np.linspace(0.026,0.030, 5)
dip1 = np.linspace(0.0005, 0.0015, 6)
dip2 = np.linspace(0.0055, 0.0075, 11)
#dip1 = np.linspace(-0.001, 0.001, 3)
#dip2 = np.linspace(0.005, 0.007, 3)

ls_div = []
ls_dip = []
ls_dth = []
p_vals_DS8 = []
z = 0
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
                # print("%d;%d;%d;%d" % (d1,d2,p1,p2))

dict = {'DIP rate': ls_dip, 'division rate': ls_div,
        'death rate': ls_dth, 'p-value': p_vals_DS8}

df = pd.DataFrame(data=dict)

df.to_pickle('PC9-DS8_param-scan_twoState_lowMedDivsMoreDIP.pkl')
