#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script of a simple 3 box model and plotting output (from Koven et al. 2015).
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

from pylab import title,savefig,figure,mean,arange
from pylab import subplot,subplots_adjust,xlim,ylim,axhline


#%%

# Length and timestep of run
nyrs = 70
timestep = 1.0
time = timestep*arange(nyrs)

# Parameters of box model
e_box = [0.3, 0.3, 0.0]
tau_box = [1., 10., 100.]
mbox = len(tau_box)

# Initial NPP
NPP = np.zeros(nyrs)
NPP[0] = 50.0

# Percentage increase in NPP per year
pc_per_yr = 0.3
mult = 1 + 0.01*pc_per_yr

# Other Array definitions
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


for n in range(1, nyrs):
    NPP[n] = NPP[n - 1]*mult 
    dCdt_box[0, n] = NPP[n] - C_box[0, n-1]/tau_box[0]
    C_box[0, n] = C_box[0, n-1] + dCdt_box[0, n]*timestep
    Rh[n] = (1 - e_box[0])*C_box[0, n]/tau_box[0]
    for m in range(1, mbox):
        dCdt_box[m, n] = e_box[m - 1]*C_box[m - 1, n]/tau_box[m - 1] - C_box[m, n-1]/tau_box[m] 
        C_box[m, n] = C_box[m, n-1] + dCdt_box[m, n]*timestep
        Rh[n] = Rh[n] + (1 - e_box[m])*C_box[m, n]/tau_box[m]
    pass
    C_tot[n] = sum(C_box[:, n])
pass
tau = C_tot/Rh


#%%

# Length and timestep of run
nyrs2 = 500
timestep2 = 1.0
time2 = timestep2*arange(nyrs2)

# Parameters of box model
e_box2 = [0.3, 0.3, 0.0]
tau_box2 = [1., 10., 100.]
mbox2 = len(tau_box2)

# Initial NPP
NPP2 = np.zeros(nyrs2)
NPP2[0] = 50.0

# Abrupt change in NPP
NPP_max2 = 60.0
time_abrupt2 = 100.0

# Other Array definitions
Rh2 = np.zeros(nyrs2)
C_tot2 = np.zeros(nyrs2)
C_box2 = np.zeros((mbox2, nyrs2))
dCdt_box2 = np.zeros((mbox2, nyrs2))


# Initialise at equilibrium
C_box2[0, 0] = NPP2[0]*tau_box2[0]
Rh2[0] = (1 - e_box2[0])*C_box2[0, 0]/tau_box2[0]
for m in range(1, mbox2):
    C_box2[m, 0] = e_box2[m - 1]*C_box2[m - 1, 0]*tau_box2[m]/tau_box2[m - 1]
    Rh2[0] = Rh2[0] + (1 - e_box2[m])*C_box2[m, 0]/tau_box2[m]
pass
C_tot2[0] = sum(C_box2[:, 0])


for n in range(1, nyrs2):
    NPP2[n] = NPP_max2
    if time2[n] < time_abrupt2:
        NPP2[n] = NPP2[0]
    dCdt_box2[0, n] = NPP2[n] - C_box2[0, n-1]/tau_box2[0]
    C_box2[0, n] = C_box2[0, n-1] + dCdt_box2[0, n]*timestep2
    Rh2[n] = (1 - e_box2[0])*C_box2[0, n]/tau_box2[0]
    for m in range(1, mbox2):
        dCdt_box2[m, n] = e_box2[m - 1]*C_box2[m - 1, n]/tau_box2[m - 1] - C_box2[m, n - 1]/tau_box2[m] 
        C_box2[m, n] = C_box2[m, n-1] + dCdt_box2[m, n]*timestep2
        Rh2[n] = Rh2[n] + (1 - e_box2[m])*C_box2[m, n]/tau_box2[m]
    pass
    C_tot2[n] = sum(C_box2[:, n])
pass
tau2 = C_tot2/Rh2


# Diagnose dCs,npp and dCs,tau
dCsnpp2 = tau2[0]*(NPP2[:] - NPP2[0])
dCstau2 = NPP2[0]*(tau2[:] - tau2[0])
print('dCstau/dCsnpp (at end) = ', dCstau2[nyrs2 - 1]/dCsnpp2[nyrs2 - 1])


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(32,26))
gs = gspec.GridSpec(2, 2, figure=fig_figure1, hspace=0.3, wspace=0.3)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':8,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':48,
    'xtick.labelsize':48,
    'ytick.labelsize':48,
    'font.size':48,
}
plt.rcParams.update(params)


## Plotting

ax = fig_figure1.add_subplot(gs[0, 0])

plt.ylabel('Normalised values')
plt.xlabel('Time (yr)')
plt.title('(a)', x=-0.3, y=1.07, fontweight='bold')

y = tau/tau[0]
y1 = NPP/NPP[0]
y2 = Rh/Rh[0]
y3 = C_tot/C_tot[0]
ylim(0.85, 1.25)
plt.plot(time, y1, 'g-', label='NPP')
plt.plot(time, y2, '-', color='lightsalmon', label='$R_h$')
plt.plot(time, y3, 'b-.', label='Total C$_s$')
plt.plot(time, y, 'r-', label=r'$\tau_s$')
plt.axhline(1, color='black', linewidth=0.5)


ax = fig_figure1.add_subplot(gs[0, 1])
colors = ['powderblue', 'mediumturquoise', 'dodgerblue']
plt.ylabel(r'Fractional $\Delta C_{s}$')
plt.xlabel('Time (yr)')

for m in range(0, mbox):
    lab = 'Box '+str(m+1)
    y = C_box[m, :]/C_box[m, 0] - 1
    plt.plot(time, y, color=colors[m])
pass
y = C_tot/C_tot[0] - 1
plt.plot(time,y,'b-.')


ax = fig_figure1.add_subplot(gs[1, 0])

plt.ylabel('Normalised values')
plt.xlabel('Time (yr)')
plt.title('(b)', x=-0.3, y=1.07, fontweight='bold')

y_2 = tau2/tau2[0]
y1_2 = NPP2/NPP2[0]
y2_2 = Rh2/Rh2[0]
y3_2 = C_tot2/C_tot2[0]
ylim(0.85, 1.25)
plt.plot(time2, y1_2, 'g-', label='NPP')
plt.plot(time2, y2_2, '-', color='lightsalmon', label='$R_h$')
plt.plot(time2, y3_2, 'b-.', label='Total C$_s$')
plt.plot(time2, y_2, 'r-', label=r'$\tau_s$')
plt.axhline(1, color='black', linewidth=0.5)


handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='g', label='NPP')])
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='lightsalmon', label=r'$R_{h}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='r', label=r'$\tau_{s}$')])
handles.extend([Line2D([0,0],[0,0], linestyle='-.', linewidth=8, color='b', label=r'Total C$_{s}$')])
labels = ['NPP', r'$R_{h}$', r'$\tau_{s}$', r'Total $C_{s}$',]
leg1 = ax.legend(handles, labels, loc='lower center', ncol=2, borderaxespad=0.2, bbox_to_anchor=(0.5, -0.6), fontsize=48)
plt.gca().add_artist(leg1)


ax = fig_figure1.add_subplot(gs[1, 1])
colors = ['powderblue', 'mediumturquoise', 'dodgerblue']
plt.ylabel(r'Fractional $\Delta C_{s}$')
plt.xlabel('Time (yr)')

for m in range(0, mbox2):
    lab = 'Box '+str(m+1)
    y_2 = C_box2[m, :]/C_box2[m, 0] - 1
    plt.plot(time2, y_2, color=colors[m])
pass
y_2 = C_tot2/C_tot2[0] - 1
plt.plot(time2, y_2, 'b-.')


handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='powderblue', label='$C_{s, 1}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='mediumturquoise', label='$C_{s, 2}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=8, color='dodgerblue', label='$C_{s, 3}$')])
handles.extend([Line2D([0,0],[0,0], linestyle='-.', linewidth=8, color='b', label='Total $C_{s}$')])
labels = ['$C_{s, 1}$', '$C_{s, 2}$', '$C_{s, 3}$', 'Total $C_{s}$']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=2, borderaxespad=0.2, bbox_to_anchor=(0.5, -0.6), fontsize=48)
plt.gca().add_artist(leg1)


#%%
fig_figure1.savefig('figures/Koven2015_subplot_v1', bbox_inches='tight')
plt.close()