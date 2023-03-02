#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots barchart of dCs and the breakdown values (dCs,npp / dCs,tau / dnppdtau, dCs,nep / dCs,tau_nep / dnepdtau)
for the C4MIP simulations for 2xCO2 and 4xCO2.
"

"""

#%%

# Analysis imports
import numpy as np
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D
#import matplotlib
#matplotlib.use('Agg')


#%% Loading in correlation arrays

# 1pctCO2 (2xCO2)
cmip6_deltaCs_1pctCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_2xCO2.npy')
cmip6_deltaCstau_1pctCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_2xCO2.npy')
cmip6_deltadelta_1pctCO2 = np.load('saved_data/cmip6_deltadelta_1pctCO2_2xCO2.npy')
cmip6_NEPtau_1pctCO2 = np.load('saved_data/cmip6_NEPtau_1pctCO2_2xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2 = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_2xCO2.npy')
cmip6_NEPdeltatau_1pctCO2 = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_2xCO2.npy')
# 1pctCO2-bgc (2xCO2)
cmip6_deltaCs_1pctCO2_bgc = np.load('saved_data/cmip6_Cs_1pctCO2_bgc_2xCO2.npy')
cmip6_deltaCstau_1pctCO2_bgc = np.load('saved_data/cmip6_Cstau_1pctCO2_bgc_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_bgc = np.load('saved_data/cmip6_CsNPP_1pctCO2_bgc_2xCO2.npy')
cmip6_deltadelta_1pctCO2_bgc = np.load('saved_data/cmip6_deltadelta_1pctCO2_bgc_2xCO2.npy')
cmip6_NEPtau_1pctCO2_bgc = np.load('saved_data/cmip6_NEPtau_1pctCO2_bgc_2xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2_bgc = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_bgc_2xCO2.npy')
cmip6_NEPdeltatau_1pctCO2_bgc = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_bgc_2xCO2.npy')
# 1pctCO2-rad (2xCO2)
cmip6_deltaCs_1pctCO2_rad = np.load('saved_data/cmip6_Cs_1pctCO2_rad_2xCO2.npy')
cmip6_deltaCstau_1pctCO2_rad = np.load('saved_data/cmip6_Cstau_1pctCO2_rad_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_rad = np.load('saved_data/cmip6_CsNPP_1pctCO2_rad_2xCO2.npy')
cmip6_deltadelta_1pctCO2_rad = np.load('saved_data/cmip6_deltadelta_1pctCO2_rad_2xCO2.npy')
cmip6_NEPtau_1pctCO2_rad = np.load('saved_data/cmip6_NEPtau_1pctCO2_rad_2xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2_rad = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_rad_2xCO2.npy')
cmip6_NEPdeltatau_1pctCO2_rad = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_rad_2xCO2.npy')


# 1pctCO2 (4xCO2)
cmip6_deltaCs_1pctCO2_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_4xCO2.npy')
cmip6_deltadelta_1pctCO2_4xCO2 = np.load('saved_data/cmip6_deltadelta_1pctCO2_4xCO2.npy')
cmip6_NEPtau_1pctCO2_4xCO2 = np.load('saved_data/cmip6_NEPtau_1pctCO2_4xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2_4xCO2 = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_4xCO2.npy')
cmip6_NEPdeltatau_1pctCO2_4xCO2 = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_4xCO2.npy')
# 1pctCO2-bgc (4xCO2)
cmip6_deltaCs_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_bgc_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_bgc_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_bgc_4xCO2.npy')
cmip6_deltadelta_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_deltadelta_1pctCO2_bgc_4xCO2.npy')
cmip6_NEPtau_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_NEPtau_1pctCO2_bgc_4xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_bgc_4xCO2.npy')
cmip6_NEPdeltatau_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_bgc_4xCO2.npy')
# 1pctCO2-rad (4xCO2)
cmip6_deltaCs_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_rad_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_rad_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_rad_4xCO2.npy')
cmip6_deltadelta_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_deltadelta_1pctCO2_rad_4xCO2.npy')
cmip6_NEPtau_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_NEPtau_1pctCO2_rad_4xCO2.npy')
cmip6_deltaNEPdeltatau_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_deltaNEPdeltatau_1pctCO2_rad_4xCO2.npy')
cmip6_NEPdeltatau_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_NEPdeltatau_1pctCO2_rad_4xCO2.npy')


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(62,96))
gs = gspec.GridSpec(3, 2, figure=fig_figure1, hspace=1.15, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':5,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':84,
    'xtick.labelsize':84,
    'ytick.labelsize':84,
    'font.size':84,
}
plt.rcParams.update(params)

#
label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']


#%% 2xCO2

#%%
column = 0
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars


cmip6_NEPtau_1pctCO2 = np.negative(cmip6_NEPtau_1pctCO2)
cmip6_deltaNEPdeltatau_1pctCO2 = np.negative(cmip6_deltaNEPdeltatau_1pctCO2)
cmip6_NEPdeltatau_1pctCO2 = np.negative(cmip6_NEPdeltatau_1pctCO2)

arrays = [cmip6_deltaCsNPP_1pctCO2, cmip6_deltaCstau_1pctCO2, cmip6_deltadelta_1pctCO2, cmip6_NEPtau_1pctCO2, cmip6_NEPdeltatau_1pctCO2, cmip6_deltaNEPdeltatau_1pctCO2]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]
    
#
ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-1750, 1750))
ax.set_title(r'(a) 2xCO$_{2}$', y=1.15, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% legend
handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label=r'$\Delta C_{s}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='g', label=r'$\Delta C_{s, NPP}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='r', label=r'$\Delta C_{s, \tau}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='k', label=r'delta NPP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label=r'\Delta NEP tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='pink', label=r'NEP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='grey', label=r'delta NEP delta tau')])
labels = [r'$\Delta C_{s}$', r'$\Delta C_{s, NPP}$', r'$\Delta C_{s, \tau}$', r'$\Delta$ NPP $\Delta \tau$', r'$\Delta C_{s, NEP}$', r'$\Delta C_{s, \tau_{NEP}}$', r'$\Delta$ NEP $\Delta \tau$']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=7, bbox_to_anchor=(1, -0.95), fontsize=84)
plt.gca().add_artist(leg1)


#%%
column = 0
row = 1
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_1pctCO2_bgc = np.negative(cmip6_NEPtau_1pctCO2_bgc)
cmip6_deltaNEPdeltatau_1pctCO2_bgc = np.negative(cmip6_deltaNEPdeltatau_1pctCO2_bgc)
cmip6_NEPdeltatau_1pctCO2_bgc = np.negative(cmip6_NEPdeltatau_1pctCO2_bgc)

arrays = [cmip6_deltaCsNPP_1pctCO2_bgc, cmip6_deltaCstau_1pctCO2_bgc, cmip6_deltadelta_1pctCO2_bgc, cmip6_NEPtau_1pctCO2_bgc, cmip6_NEPdeltatau_1pctCO2_bgc, cmip6_deltaNEPdeltatau_1pctCO2_bgc]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_bgc, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2_bgc))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2_bgc))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}^{BGC}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-1750, 1750))
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% legend
handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label=r'$\Delta C_{s}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='g', label=r'$\Delta C_{s, NPP}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='r', label=r'$\Delta C_{s, \tau}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='k', label=r'delta NPP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label=r'\Delta NEP tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='pink', label=r'NEP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='grey', label=r'delta NEP delta tau')])
labels = [r'$\Delta C_{s}^{BGC}$', r'$\Delta C_{s, NPP}^{BGC}$', r'$\Delta C_{s, \tau}^{BGC}$', r'$\Delta NPP^{BGC}$ $\Delta \tau^{BGC}$', r'$\Delta C_{s, NEP}^{BGC}$', r'$\Delta C_{s, \tau_{NEP}}^{BGC}$', r'$\Delta NEP^{BGC}$ $\Delta \tau^{BGC}$']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=7, bbox_to_anchor=(1, -0.95), fontsize=84)
plt.gca().add_artist(leg1)


#%%
column = 0
row = 2
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_1pctCO2_rad = np.negative(cmip6_NEPtau_1pctCO2_rad)
cmip6_deltaNEPdeltatau_1pctCO2_rad = np.negative(cmip6_deltaNEPdeltatau_1pctCO2_rad)
cmip6_NEPdeltatau_1pctCO2_rad = np.negative(cmip6_NEPdeltatau_1pctCO2_rad)

arrays = [cmip6_deltaCsNPP_1pctCO2_rad, cmip6_deltaCstau_1pctCO2_rad, cmip6_deltadelta_1pctCO2_rad, cmip6_NEPtau_1pctCO2_rad, cmip6_NEPdeltatau_1pctCO2_rad, cmip6_deltaNEPdeltatau_1pctCO2_rad]
colours = ['g', 'r',  'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_rad, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2_rad))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2_rad))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}^{RAD}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-250, 250))
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% legend
handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label=r'$\Delta C_{s}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='g', label=r'$\Delta C_{s, NPP}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='r', label=r'$\Delta C_{s, \tau}$')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='k', label=r'delta NPP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label=r'\Delta NEP tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='pink', label=r'NEP delta tau')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='grey', label=r'delta NEP delta tau')])
labels = [r'$\Delta C_{s}^{RAD}$', r'$\Delta C_{s, NPP}^{RAD}$', r'$\Delta C_{s, \tau}^{RAD}$', r'$\Delta NPP^{RAD}$ $\Delta \tau^{RAD}$', r'$\Delta C_{s, NEP}^{RAD}$', r'$\Delta C_{s, \tau_{NEP}}^{RAD}$', r'$\Delta NEP^{RAD}$ $\Delta \tau^{RAD}$']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=7, bbox_to_anchor=(1, -0.95), fontsize=84)#title='Model colors',
plt.gca().add_artist(leg1)


#%% 4xCO2 

#%%
column = 1
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars


cmip6_NEPtau_1pctCO2_4xCO2 = np.negative(cmip6_NEPtau_1pctCO2_4xCO2)
cmip6_deltaNEPdeltatau_1pctCO2_4xCO2 = np.negative(cmip6_deltaNEPdeltatau_1pctCO2_4xCO2)
cmip6_NEPdeltatau_1pctCO2_4xCO2 = np.negative(cmip6_NEPdeltatau_1pctCO2_4xCO2)

arrays = [cmip6_deltaCsNPP_1pctCO2_4xCO2, cmip6_deltaCstau_1pctCO2_4xCO2, cmip6_deltadelta_1pctCO2_4xCO2, cmip6_NEPtau_1pctCO2_4xCO2, cmip6_NEPdeltatau_1pctCO2_4xCO2, cmip6_deltaNEPdeltatau_1pctCO2_4xCO2]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_4xCO2, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2_4xCO2))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2_4xCO2))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]
    
#
ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-3250, 3250))
ax.set_title(r'(b) 4xCO$_{2}$', y=1.15, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%%
column = 1
row = 1
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_1pctCO2_bgc_4xCO2 = np.negative(cmip6_NEPtau_1pctCO2_bgc_4xCO2)
cmip6_deltaNEPdeltatau_1pctCO2_bgc_4xCO2 = np.negative(cmip6_deltaNEPdeltatau_1pctCO2_bgc_4xCO2)
cmip6_NEPdeltatau_1pctCO2_bgc_4xCO2 = np.negative(cmip6_NEPdeltatau_1pctCO2_bgc_4xCO2)

arrays = [cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2, cmip6_deltaCstau_1pctCO2_bgc_4xCO2, cmip6_deltadelta_1pctCO2_bgc_4xCO2, cmip6_NEPtau_1pctCO2_bgc_4xCO2, cmip6_NEPdeltatau_1pctCO2_bgc_4xCO2, cmip6_deltaNEPdeltatau_1pctCO2_bgc_4xCO2]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_bgc_4xCO2, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2_bgc_4xCO2))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2_bgc_4xCO2))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}^{BGC}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-3250, 3250))
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%%
column = 1
row = 2
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_1pctCO2_rad_4xCO2 = np.negative(cmip6_NEPtau_1pctCO2_rad_4xCO2)
cmip6_deltaNEPdeltatau_1pctCO2_rad_4xCO2 = np.negative(cmip6_deltaNEPdeltatau_1pctCO2_rad_4xCO2)
cmip6_NEPdeltatau_1pctCO2_rad_4xCO2 = np.negative(cmip6_NEPdeltatau_1pctCO2_rad_4xCO2)

arrays = [cmip6_deltaCsNPP_1pctCO2_rad_4xCO2, cmip6_deltaCstau_1pctCO2_rad_4xCO2, cmip6_deltadelta_1pctCO2_rad_4xCO2, cmip6_NEPtau_1pctCO2_rad_4xCO2, cmip6_NEPdeltatau_1pctCO2_rad_4xCO2, cmip6_deltaNEPdeltatau_1pctCO2_rad_4xCO2]
colours = ['g', 'r',  'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_rad_4xCO2, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_1pctCO2_rad_4xCO2))
bottom_pos = np.zeros(len(cmip6_deltaCs_1pctCO2_rad_4xCO2))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}^{RAD}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-650, 650))
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%%
fig_figure1.savefig('figures/dCs_breakdown_barchart_C4MIP_v2', bbox_inches='tight')
plt.close()