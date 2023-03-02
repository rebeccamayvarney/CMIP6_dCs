#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots barchart of dCs and the breakdown values (dCs,npp / dCs,tau / dnppdtau, dCs,nep / dCs,tau_nep / dnepdtau)
for CMIP5 and CMIP6 ESMs for RCP and SSP simulations.
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


#%% Loading in correlation arrays

# SSP126 (CMIP6)
cmip6_deltaCs_ssp126 = np.load('saved_data/cmip6_Cs_ssp126_npp.npy')
cmip6_deltaCstau_ssp126 = np.load('saved_data/cmip6_Cstau_ssp126_npp.npy')
cmip6_deltaCsRh_ssp126 = np.load('saved_data/cmip6_CsRh_ssp126_npp.npy')
cmip6_deltadelta_ssp126 = np.load('saved_data/cmip6_deltadelta_ssp126_npp.npy')
cmip6_remaining_ssp126 = np.load('saved_data/cmip6_remaining_ssp126_npp.npy')
cmip6_NEPtau_ssp126 = np.load('saved_data/cmip6_NEPtau_ssp126.npy')
cmip6_deltaNEPdeltatau_ssp126 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp126.npy')
cmip6_NEPdeltatau_ssp126 = np.load('saved_data/cmip6_NEPdeltatau_ssp126.npy')
# RCP2.6 (CMIP5)
cmip5_deltaCs_rcp26 = np.load('saved_data/cmip5_Cs_rcp26_npp.npy')
cmip5_deltaCstau_rcp26 = np.load('saved_data/cmip5_Cstau_rcp26_npp.npy')
cmip5_deltaCsRh_rcp26 = np.load('saved_data/cmip5_CsRh_rcp26_npp.npy')
cmip5_deltadelta_rcp26 = np.load('saved_data/cmip5_deltadelta_rcp26_npp.npy')
cmip5_remaining_rcp26 = np.load('saved_data/cmip5_remaining_rcp26_npp.npy')
cmip5_NEPtau_rcp26 = np.load('saved_data/cmip5_NEPtau_rcp26_npp.npy')
cmip5_deltaNEPdeltatau_rcp26 = np.load('saved_data/cmip5_deltaNEPdeltatau_rcp26.npy')
cmip5_NEPdeltatau_rcp26 = np.load('saved_data/cmip5_NEPdeltatau_rcp26.npy')

# SSP245 (CMIP6)
cmip6_deltaCs_ssp245 = np.load('saved_data/cmip6_Cs_ssp245_npp.npy')
cmip6_deltaCstau_ssp245 = np.load('saved_data/cmip6_Cstau_ssp245_npp.npy')
cmip6_deltaCsRh_ssp245 = np.load('saved_data/cmip6_CsRh_ssp245_npp.npy')
cmip6_deltadelta_ssp245 = np.load('saved_data/cmip6_deltadelta_ssp245_npp.npy')
cmip6_remaining_ssp245 = np.load('saved_data/cmip6_remaining_ssp245_npp.npy')
cmip6_NEPtau_ssp245 = np.load('saved_data/cmip6_NEPtau_ssp245.npy')
cmip6_deltaNEPdeltatau_ssp245 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp245.npy')
cmip6_NEPdeltatau_ssp245 = np.load('saved_data/cmip6_NEPdeltatau_ssp245.npy')
# RCP4.5 (CMIP5)
cmip5_deltaCs_rcp45 = np.load('saved_data/cmip5_Cs_rcp45_npp.npy')
cmip5_deltaCstau_rcp45 = np.load('saved_data/cmip5_Cstau_rcp45_npp.npy')
cmip5_deltaCsRh_rcp45 = np.load('saved_data/cmip5_CsRh_rcp45_npp.npy')
cmip5_deltadelta_rcp45 = np.load('saved_data/cmip5_deltadelta_rcp45_npp.npy')
cmip5_remaining_rcp45 = np.load('saved_data/cmip5_remaining_rcp45_npp.npy')
cmip5_NEPtau_rcp45 = np.load('saved_data/cmip5_NEPtau_rcp45_npp.npy')
cmip5_deltaNEPdeltatau_rcp45 = np.load('saved_data/cmip5_deltaNEPdeltatau_rcp45.npy')
cmip5_NEPdeltatau_rcp45 = np.load('saved_data/cmip5_NEPdeltatau_rcp45.npy')

# SSP585 (CMIP6)
cmip6_deltaCs_ssp585 = np.load('saved_data/cmip6_Cs_ssp585_npp.npy')
cmip6_deltaCstau_ssp585 = np.load('saved_data/cmip6_Cstau_ssp585_npp.npy')
cmip6_deltaCsRh_ssp585 = np.load('saved_data/cmip6_CsRh_ssp585_npp.npy')
cmip6_deltadelta_ssp585 = np.load('saved_data/cmip6_deltadelta_ssp585_npp.npy')
cmip6_remaining_ssp585 = np.load('saved_data/cmip6_remaining_ssp585_npp.npy')
cmip6_NEPtau_ssp585 = np.load('saved_data/cmip6_NEPtau_ssp585.npy')
cmip6_deltaNEPdeltatau_ssp585 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp585.npy')
cmip6_NEPdeltatau_ssp585 = np.load('saved_data/cmip6_NEPdeltatau_ssp585.npy')
# RCP8.5 (CMIP5)
cmip5_deltaCs_rcp85 = np.load('saved_data/cmip5_Cs_rcp85_npp.npy')
cmip5_deltaCstau_rcp85 = np.load('saved_data/cmip5_Cstau_rcp85_npp.npy')
cmip5_deltaCsRh_rcp85 = np.load('saved_data/cmip5_CsRh_rcp85_npp.npy')
cmip5_deltadelta_rcp85 = np.load('saved_data/cmip5_deltadelta_rcp85_npp.npy')
cmip5_remaining_rcp85 = np.load('saved_data/cmip5_remaining_rcp85_npp.npy')
cmip5_NEPtau_rcp85 = np.load('saved_data/cmip5_NEPtau_rcp85_npp.npy')
cmip5_deltaNEPdeltatau_rcp85 = np.load('saved_data/cmip5_deltaNEPdeltatau_rcp85.npy')
cmip5_NEPdeltatau_rcp85 = np.load('saved_data/cmip5_NEPdeltatau_rcp85.npy')


#%% Merging CMIP5 & CMIP6

# dCs
cmip5_Cs = np.concatenate([cmip5_deltaCs_rcp26, cmip5_deltaCs_rcp45, cmip5_deltaCs_rcp85])
cmip6_Cs = np.concatenate([cmip6_deltaCs_ssp126, cmip6_deltaCs_ssp245, cmip6_deltaCs_ssp585])
# dCs,npp
cmip5_CsNPP = np.concatenate([cmip5_deltaCsRh_rcp26, cmip5_deltaCsRh_rcp45, cmip5_deltaCsRh_rcp85])
cmip6_CsNPP = np.concatenate([cmip6_deltaCsRh_ssp126, cmip6_deltaCsRh_ssp245, cmip6_deltaCsRh_ssp585])
# dCs,tau
cmip5_CsTAU = np.concatenate([cmip5_deltaCstau_rcp26, cmip5_deltaCstau_rcp45, cmip5_deltaCstau_rcp85])
cmip6_CsTAU = np.concatenate([cmip6_deltaCstau_ssp126, cmip6_deltaCstau_ssp245, cmip6_deltaCstau_ssp585])


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(48,38))
gs = gspec.GridSpec(2, 3, figure=fig_figure1, hspace=0.5, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':52,
    'xtick.labelsize':52,
    'ytick.labelsize':52,
    'font.size':52,
}
plt.rcParams.update(params)


#%% CMIP6
label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']


#%% cSoil
column = 0
row = 1
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars


cmip6_NEPtau_ssp126 = np.negative(cmip6_NEPtau_ssp126)
cmip6_deltaNEPdeltatau_ssp126 = np.negative(cmip6_deltaNEPdeltatau_ssp126)
cmip6_NEPdeltatau_ssp126 = np.negative(cmip6_NEPdeltatau_ssp126)

arrays = [cmip6_deltaCsRh_ssp126, cmip6_deltaCstau_ssp126, cmip6_deltadelta_ssp126, cmip6_NEPtau_ssp126, cmip6_NEPdeltatau_ssp126, cmip6_deltaNEPdeltatau_ssp126]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_ssp126, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_ssp126))
bottom_pos = np.zeros(len(cmip6_deltaCs_ssp126))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]
    
#
ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-800, 800))
ax.set_title('SSP126', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% NPP
column = 1
row = 1
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_ssp245 = np.negative(cmip6_NEPtau_ssp245)
cmip6_deltaNEPdeltatau_ssp245 = np.negative(cmip6_deltaNEPdeltatau_ssp245)
cmip6_NEPdeltatau_ssp245 = np.negative(cmip6_NEPdeltatau_ssp245)

arrays = [cmip6_deltaCsRh_ssp245, cmip6_deltaCstau_ssp245, cmip6_deltadelta_ssp245, cmip6_NEPtau_ssp245, cmip6_NEPdeltatau_ssp245, cmip6_deltaNEPdeltatau_ssp245]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_ssp245, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_ssp245))
bottom_pos = np.zeros(len(cmip6_deltaCs_ssp245))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]

#
ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-1400, 1400))
ax.set_title('SSP245', fontweight='bold')
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
leg1 = ax.legend(handles, labels, loc='lower center', ncol=7, bbox_to_anchor=(0.4, -0.7), fontsize=52)
plt.gca().add_artist(leg1)


#%% tau
column = 2
row = 1
ax = fig_figure1.add_subplot(gs[row, column])


cmip6_NEPtau_ssp585 = np.negative(cmip6_NEPtau_ssp585)
cmip6_deltaNEPdeltatau_ssp585 = np.negative(cmip6_deltaNEPdeltatau_ssp585)
cmip6_NEPdeltatau_ssp585 = np.negative(cmip6_NEPdeltatau_ssp585)

arrays = [cmip6_deltaCsRh_ssp585, cmip6_deltaCstau_ssp585, cmip6_deltadelta_ssp585, cmip6_NEPtau_ssp585, cmip6_NEPdeltatau_ssp585, cmip6_deltaNEPdeltatau_ssp585]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_ssp585, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip6_deltaCs_ssp585))
bottom_pos = np.zeros(len(cmip6_deltaCs_ssp585))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]

#
ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-3300, 3300))
ax.set_title('SSP585', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% CMIP5
label_list = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']


#%% cSoil
column = 0
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(9)  # the label locations
width = 0.4  # the width of the bars


cmip5_NEPtau_rcp26 = np.negative(cmip5_NEPtau_rcp26)
cmip5_deltaNEPdeltatau_rcp26 = np.negative(cmip5_deltaNEPdeltatau_rcp26)
cmip5_NEPdeltatau_rcp26 = np.negative(cmip5_NEPdeltatau_rcp26)

arrays = [cmip5_deltaCsRh_rcp26, cmip5_deltaCstau_rcp26, cmip5_deltadelta_rcp26, cmip5_NEPtau_rcp26, cmip5_NEPdeltatau_rcp26, cmip5_deltaNEPdeltatau_rcp26]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip5_deltaCs_rcp26, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip5_deltaCs_rcp26))
bottom_pos = np.zeros(len(cmip5_deltaCs_rcp26))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-800, 800))
ax.set_title('RCP2.6', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%% NPP
column = 1
row = 0
ax = fig_figure1.add_subplot(gs[row, column])


cmip5_NEPtau_rcp45 = np.negative(cmip5_NEPtau_rcp45)
cmip5_deltaNEPdeltatau_rcp45 = np.negative(cmip5_deltaNEPdeltatau_rcp45)
cmip5_NEPdeltatau_rcp45 = np.negative(cmip5_NEPdeltatau_rcp45)

arrays = [cmip5_deltaCsRh_rcp45, cmip5_deltaCstau_rcp45, cmip5_deltadelta_rcp45, cmip5_NEPtau_rcp45, cmip5_NEPdeltatau_rcp45, cmip5_deltaNEPdeltatau_rcp45]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip5_deltaCs_rcp45, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip5_deltaCs_rcp45))
bottom_pos = np.zeros(len(cmip5_deltaCs_rcp45))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-1400, 1400))
ax.set_title('RCP4.5', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')



#%% tau
column = 2
row = 0
ax = fig_figure1.add_subplot(gs[row, column])


cmip5_NEPtau_rcp85 = np.negative(cmip5_NEPtau_rcp85)
cmip5_deltaNEPdeltatau_rcp85 = np.negative(cmip5_deltaNEPdeltatau_rcp85)
cmip5_NEPdeltatau_rcp85 = np.negative(cmip5_NEPdeltatau_rcp85)

arrays = [cmip5_deltaCsRh_rcp85, cmip5_deltaCstau_rcp85, cmip5_deltadelta_rcp85, cmip5_NEPtau_rcp85, cmip5_NEPdeltatau_rcp85, cmip5_deltaNEPdeltatau_rcp85]
colours = ['g', 'r', 'k', 'darkseagreen', 'pink', 'grey']

rects1 = ax.bar(x - 1/5, cmip5_deltaCs_rcp85, width, color='b', label='Cs')

bottom_neg = np.zeros(len(cmip5_deltaCs_rcp85))
bottom_pos = np.zeros(len(cmip5_deltaCs_rcp85))

for i in range(0,6):
    bottom = (arrays[i]<0)*bottom_neg + (arrays[i]>0)*bottom_pos
    ax.bar(x + 1/5, arrays[i], width, bottom = bottom, color=colours[i])
    bottom_neg = bottom_neg + (arrays[i]<0)*arrays[i]
    bottom_pos = bottom_pos + (arrays[i]>0)*arrays[i]


ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-3300, 3300))
ax.set_title('RCP8.5', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%%
fig_figure1.savefig('figures/dCs_breakdown_barchart_cmip5cmip6_v4', bbox_inches='tight')
plt.close()