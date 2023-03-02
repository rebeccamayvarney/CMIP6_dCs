#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 2023
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dCs,npp Vs dCs,tau for CMIP6 ESMs using the future SSP and RCP simulations.
Key figure for Biogeosciences submission.
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


#%% Loading in arrays

# SSP126 (CMIP6)
cmip6_deltaCs_ssp126 = np.load('saved_data/cmip6_Cs_ssp126_npp.npy')
cmip6_deltaCstau_ssp126 = np.load('saved_data/cmip6_Cstau_ssp126_npp.npy')
cmip6_deltaCsRh_ssp126 = np.load('saved_data/cmip6_CsRh_ssp126_npp.npy')
cmip6_deltadelta_ssp126 = np.load('saved_data/cmip6_deltadelta_ssp126_npp.npy')
cmip6_remaining_ssp126 = np.load('saved_data/cmip6_remaining_ssp126_npp.npy')
cmip6_NEPtau_ssp126 = np.load('saved_data/cmip6_NEPtau_ssp126.npy')
cmip6_deltaNEPdeltatau_ssp126 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp126.npy')
cmip6_NEPdeltatau_ssp126 = np.load('saved_data/cmip6_NEPdeltatau_ssp126.npy')

# SSP245 (CMIP6)
cmip6_deltaCs_ssp245 = np.load('saved_data/cmip6_Cs_ssp245_npp.npy')
cmip6_deltaCstau_ssp245 = np.load('saved_data/cmip6_Cstau_ssp245_npp.npy')
cmip6_deltaCsRh_ssp245 = np.load('saved_data/cmip6_CsRh_ssp245_npp.npy')
cmip6_deltadelta_ssp245 = np.load('saved_data/cmip6_deltadelta_ssp245_npp.npy')
cmip6_remaining_ssp245 = np.load('saved_data/cmip6_remaining_ssp245_npp.npy')
cmip6_NEPtau_ssp245 = np.load('saved_data/cmip6_NEPtau_ssp245.npy')
cmip6_deltaNEPdeltatau_ssp245 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp245.npy')
cmip6_NEPdeltatau_ssp245 = np.load('saved_data/cmip6_NEPdeltatau_ssp245.npy')

# SSP585 (CMIP6)
cmip6_deltaCs_ssp585 = np.load('saved_data/cmip6_Cs_ssp585_npp.npy')
cmip6_deltaCstau_ssp585 = np.load('saved_data/cmip6_Cstau_ssp585_npp.npy')
cmip6_deltaCsRh_ssp585 = np.load('saved_data/cmip6_CsRh_ssp585_npp.npy')
cmip6_deltadelta_ssp585 = np.load('saved_data/cmip6_deltadelta_ssp585_npp.npy')
cmip6_remaining_ssp585 = np.load('saved_data/cmip6_remaining_ssp585_npp.npy')
cmip6_NEPtau_ssp585 = np.load('saved_data/cmip6_NEPtau_ssp585.npy')
cmip6_deltaNEPdeltatau_ssp585 = np.load('saved_data/cmip6_deltaNEPdeltatau_ssp585.npy')
cmip6_NEPdeltatau_ssp585 = np.load('saved_data/cmip6_NEPdeltatau_ssp585.npy')


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(16,14))
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':5,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':32,
    'xtick.labelsize':32,
    'ytick.labelsize':32,
    'font.size':32,
}
plt.rcParams.update(params)


#%% CMIP6

# ssp scenarios
ssp_options = ['ssp126', 'ssp245', 'ssp585']
ssp_options_length = len(ssp_options)
ssp_titles = ['SSP126', 'SSP245', 'SSP585']
ssp_markers = ['o', 'v', '*']

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors_cmip6 = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'red', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


# for loop for each CMIP6 model
for model_i in range(0, n_models):
    model = cmip6_models[model_i] # seleting the models

    plt.plot(cmip6_deltaCstau_ssp126[model_i], cmip6_deltaCsRh_ssp126[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[0], markersize=30)
    plt.plot(cmip6_deltaCstau_ssp245[model_i], cmip6_deltaCsRh_ssp245[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[1], markersize=30)
    plt.plot(cmip6_deltaCstau_ssp585[model_i], cmip6_deltaCsRh_ssp585[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[2], markersize=30)
    
plt.xlabel(r'Change in soil carbon due to $\Delta \tau_{s}$')
plt.ylabel(r'Change in soil carbon due to $\Delta \rm{NPP}$')
                      
        
#%%

handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='peachpuff', label='ACCESS-ESM1-5')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#fb8072', label='BCC-CSM2-MR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM5')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='dodgerblue', label='CESM2')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='red', label='CNRM-ESM2-1')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM6A-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='olive', label='MIROC-ES2L')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='gold', label='MPI-ESM1-2-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='orange', label='NorESM2-LM')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label='UKESM1-0-LL')])
labels = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
leg1 = plt.legend(handles, labels, loc='center right', borderaxespad=0.2, bbox_to_anchor=(1.5, 0.5), title='CMIP6', fontsize=32)
plt.gca().add_artist(leg1)

handles2 = []
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=30, color='k', label='SSP126')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='v', markersize=30, color='k', label='SSP245')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=30, color='k', label='SSP585')])
labels2 = ['SSP126', 'SSP245', 'SSP585']
leg2 = plt.legend(handles2, labels2, loc='upper right', borderaxespad=0.2, fontsize=32)
plt.gca().add_artist(leg2)


#%%
fig_figure1.savefig('figures/keyfigure_v1', bbox_inches='tight')
plt.close()