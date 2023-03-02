#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dCs,npp Vs dCs,tau derived using the simple 3 box model.
"

"""


#%%

import numpy as np
# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams, colors
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D

from pylab import plot,title,savefig,figure,mean,arange
from pylab import subplot,subplots_adjust,xlim,ylim,axhline


#%%

# Length and timestep of run
nyrs = 70
timestep = 1.0
time = timestep*arange(nyrs)


# Parameters of box model
e_box = [0.3, 0.3, 0.0]
tau_box = [1., 10., 100.]
#e_box = [0.3, 0.0]
#tau_box = [1., 10000.]
mbox = len(tau_box)

# Percentage increase in NPP per year
pc_per_yr = [0.0001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8]
kmax = len(pc_per_yr)
dCsnpp = np.zeros(kmax)
dCstau = np.zeros(kmax)

# Initial NPP
NPP = np.zeros(nyrs)
NPP[0] = 50.0

# Empty arrays
Rh = np.zeros(nyrs)
C_tot = np.zeros(nyrs)
C_box = np.zeros((mbox, nyrs))
dCdt_box = np.zeros((mbox, nyrs))


# Initialise at equilibrium
C_box[0, 0] = NPP[0]*tau_box[0]
Rh[0] = (1 - e_box[0])*C_box[0, 0]/tau_box[0]

for m in range(1, mbox):
    C_box[m, 0] = e_box[m - 1]*C_box[m - 1, 0]*tau_box[m]/tau_box[m - 1]
    Rh[0] = Rh[0] + (1 - e_box[m])*C_box[m, 0]/tau_box[m]
pass
C_tot[0] = sum(C_box[:, 0])


for k in range(0, kmax):
    mult = 1 + 0.01*pc_per_yr[k]
    for n in range(1, nyrs):
        NPP[n] = NPP[n - 1]*mult 
        dCdt_box[0, n] = NPP[n] - C_box[0, n - 1]/tau_box[0]
        C_box[0, n] = C_box[0, n - 1] + dCdt_box[0, n]*timestep
        Rh[n] = (1 - e_box[0])*C_box[0, n]/tau_box[0]
        for m in range(1, mbox):
            dCdt_box[m, n] = e_box[m - 1]*C_box[m - 1, n]/tau_box[m - 1] - C_box[m, n - 1]/tau_box[m] 
            C_box[m, n] = C_box[m, n - 1] + dCdt_box[m, n]*timestep
            Rh[n] = Rh[n] + (1 - e_box[m])*C_box[m, n]/tau_box[m]
        pass
        C_tot[n] = sum(C_box[:, n])
    pass
    tau = C_tot/Rh
    
    
# Diagnose dCs,npp and dCs,tau
    dCsnpp[k] = tau[0]*(NPP[nyrs - 1] - NPP[0])
    dCstau[k] = NPP[0]*(tau[nyrs - 1] - tau[0])
    print('% increase in NPP per year = ', pc_per_yr[k], ', dCstau/dCsnpp (at end) = ', dCstau[k]/dCsnpp[k])
pass


#%% Plotting

fig_figure1 = plt.figure(1, figsize=(8,6))
#gs = gspec.GridSpec(1, 2, figure=fig_figure1, hspace=0.4, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
#    'lines.linewidth':8,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':22,
    'xtick.labelsize':22,
    'ytick.labelsize':22,
    'font.size':22,
}
plt.rcParams.update(params)

plt.xlabel(r'$\Delta C_{s, \tau}$ (PgC)')
plt.ylabel(r'$\Delta C_{s, NPP}$ (PgC)')

plt.plot(dCstau, dCsnpp, 'k-o', linewidth=3)
plt.xlim((-140,0))
plt.ylim(0,500)


fig_figure1.savefig('figures/dCsnpp_v_dCstau_simple3boxmodel_v1', bbox_inches='tight')
plt.close()
