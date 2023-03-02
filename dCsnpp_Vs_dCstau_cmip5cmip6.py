#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dCs,npp Vs dCs,tau for CMIP5 and CMIP6 ESMs using the future SSP and RCP simulations.
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


# SSP126 (CMIP6)
cmip6_fractionaldeltaCstau_ssp126 = np.load('saved_data/cmip6_fractionalCstau_ssp126_npp.npy')
cmip6_fractionaldeltaCsRh_ssp126 = np.load('saved_data/cmip6_fractionalCsRh_ssp126_npp.npy')
# RCP2.6 (CMIP5)
cmip5_fractionaldeltaCstau_rcp26 = np.load('saved_data/cmip5_fractionalCstau_rcp26_npp.npy')
cmip5_fractionaldeltaCsRh_rcp26 = np.load('saved_data/cmip5_fractionalCsRh_rcp26_npp.npy')

# SSP245 (CMIP6)
cmip6_fractionaldeltaCstau_ssp245 = np.load('saved_data/cmip6_fractionalCstau_ssp245_npp.npy')
cmip6_fractionaldeltaCsRh_ssp245 = np.load('saved_data/cmip6_fractionalCsRh_ssp245_npp.npy')
# RCP4.5 (CMIP5)
cmip5_fractionaldeltaCstau_rcp45 = np.load('saved_data/cmip5_fractionalCstau_rcp45_npp.npy')
cmip5_fractionaldeltaCsRh_rcp45 = np.load('saved_data/cmip5_fractionalCsRh_rcp45_npp.npy')

# SSP585 (CMIP6)
cmip6_fractionaldeltaCstau_ssp585 = np.load('saved_data/cmip6_fractionalCstau_ssp585_npp.npy')
cmip6_fractionaldeltaCsRh_ssp585 = np.load('saved_data/cmip6_fractionalCsRh_ssp585_npp.npy')
# RCP8.5 (CMIP5)
cmip5_fractionaldeltaCstau_rcp85 = np.load('saved_data/cmip5_fractionalCstau_rcp85_npp.npy')
cmip5_fractionaldeltaCsRh_rcp85 = np.load('saved_data/cmip5_fractionalCsRh_rcp85_npp.npy')


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(30,24))
gs = gspec.GridSpec(2, 2, figure=fig_figure1, hspace=0.25, wspace=0.25)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':5,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':40,
    'xtick.labelsize':40,
    'ytick.labelsize':40,
    'font.size':40,
}
plt.rcParams.update(params)


#%% CMIP5

# rcp scenarios
rcp_options = ['rcp26', 'rcp45', 'rcp85']
rcp_options_length = len(rcp_options)
rcp_titles = ['RCP2.6', 'RCP4.5', 'RCP8.5']
rcp_markers = ['o', 'v', '*']

# CMIP5 models
cmip5_models = cmip5_models = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)
model_colors_cmip5 = ['darkblue', '#80b1d3', 'darkcyan', '#8dd3c7', 'darkseagreen', 'darkgreen', 'olive', 'gold', 'orange']

                      
ax = fig_figure1.add_subplot(gs[0, 0])
# for loop for each CMIP5 model
for model_i in range(0, n_models):
    model = cmip5_models[model_i] # seleting the models

    plt.plot(cmip5_deltaCstau_rcp26[model_i], cmip5_deltaCsRh_rcp26[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[0], markersize=30)
    plt.plot(cmip5_deltaCstau_rcp45[model_i], cmip5_deltaCsRh_rcp45[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[1], markersize=30)
    plt.plot(cmip5_deltaCstau_rcp85[model_i], cmip5_deltaCsRh_rcp85[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[2], markersize=30)
    
plt.xlabel(r'$\Delta C_{s, \tau}$')
plt.ylabel(r'$\Delta C_{s, NPP}$')
plt.title('(a) Absolute', y=1.08, fontweight='bold')


ax = fig_figure1.add_subplot(gs[0, 1])
# for loop for each CMIP5 model
for model_i in range(0, n_models):
    model = cmip5_models[model_i] # seleting the models

    plt.plot(cmip5_fractionaldeltaCstau_rcp26[model_i], cmip5_fractionaldeltaCsRh_rcp26[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[0], markersize=30)
    plt.plot(cmip5_fractionaldeltaCstau_rcp45[model_i], cmip5_fractionaldeltaCsRh_rcp45[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[1], markersize=30)
    plt.plot(cmip5_fractionaldeltaCstau_rcp85[model_i], cmip5_fractionaldeltaCsRh_rcp85[model_i], color=model_colors_cmip5[model_i], marker=rcp_markers[2], markersize=30)
    
plt.xlabel(r'$\Delta C_{s, \tau} / C_{s,0}$')
plt.ylabel(r'$\Delta C_{s, NPP} / C_{s,0}$')
plt.title('(b) Fractional', y=1.08, fontweight='bold')

##

#cmip5_deltaCstau_rcp26 = np.delete(cmip5_deltaCstau_rcp26, 1)
#cmip5_deltaCstau_rcp45 = np.delete(cmip5_deltaCstau_rcp45, 1)
#cmip5_deltaCstau_rcp85 = np.delete(cmip5_deltaCstau_rcp85, 1)
#cmip5_deltaCsRh_rcp26 = np.delete(cmip5_deltaCsRh_rcp26, 1)
#cmip5_deltaCsRh_rcp45 = np.delete(cmip5_deltaCsRh_rcp45, 1)
#cmip5_deltaCsRh_rcp85 = np.delete(cmip5_deltaCsRh_rcp85, 1)

r_cmip5 = np.corrcoef([cmip5_deltaCstau_rcp26.flatten(), cmip5_deltaCstau_rcp45.flatten(), cmip5_deltaCstau_rcp85.flatten()], [cmip5_deltaCsRh_rcp26.flatten(), cmip5_deltaCsRh_rcp45.flatten(), cmip5_deltaCsRh_rcp85.flatten()])[0, 1]**2
print('CMIP5 delta', r_cmip5)
r_cmip5 = np.corrcoef(cmip5_fractionaldeltaCstau_rcp45.flatten(), cmip5_fractionaldeltaCsRh_rcp45.flatten())[0, 1]**2
print('CMIP5 frac', r_cmip5)     


#%%

handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkblue', label='BNU-ESM')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM2')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkcyan', label='GFDL-ESM2G')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#8dd3c7', label='GISS-E2-R')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label='HadGEM2-ES')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM5A-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='olive', label='MIROC-ESM')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='gold', label='MPI-ESM-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='orange', label='NorESM1-M')])
labels = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
leg1 = ax.legend(handles, labels, loc='center right', borderaxespad=0.2, bbox_to_anchor=(1.77, 0.5), title='CMIP5', fontsize=40)
plt.gca().add_artist(leg1)


handles2 = []
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=30, color='k', label='RCP2.6')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='v', markersize=30, color='k', label='RCP4.5')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=30, color='k', label='RCP8.5')])
labels2 = ['RCP2.6', 'RCP4.5', 'RCP8.5']
leg2 = ax.legend(handles2, labels2, loc='upper right', borderaxespad=0.2, fontsize=40)
plt.gca().add_artist(leg2)


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


ax = fig_figure1.add_subplot(gs[1, 0])
# for loop for each CMIP6 model
for model_i in range(0, n_models):
    model = cmip6_models[model_i] # seleting the models

    plt.plot(cmip6_deltaCstau_ssp126[model_i], cmip6_deltaCsRh_ssp126[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[0], markersize=30)
    plt.plot(cmip6_deltaCstau_ssp245[model_i], cmip6_deltaCsRh_ssp245[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[1], markersize=30)
    plt.plot(cmip6_deltaCstau_ssp585[model_i], cmip6_deltaCsRh_ssp585[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[2], markersize=30)
    
plt.xlabel(r'$\Delta C_{s, \tau}$')
plt.ylabel(r'$\Delta C_{s, NPP}$')
                      

ax = fig_figure1.add_subplot(gs[1, 1])
# for loop for each CMIP6 model
for model_i in range(0, n_models):
    model = cmip6_models[model_i] # seleting the models

    plt.plot(cmip6_fractionaldeltaCstau_ssp126[model_i], cmip6_fractionaldeltaCsRh_ssp126[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[0], markersize=30)
    plt.plot(cmip6_fractionaldeltaCstau_ssp245[model_i], cmip6_fractionaldeltaCsRh_ssp245[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[1], markersize=30)
    plt.plot(cmip6_fractionaldeltaCstau_ssp585[model_i], cmip6_fractionaldeltaCsRh_ssp585[model_i], color=model_colors_cmip6[model_i], marker=ssp_markers[2], markersize=30)
    
plt.xlabel(r'$\Delta C_{s, \tau} / C_{s,0}$')
plt.ylabel(r'$\Delta C_{s, NPP} / C_{s,0}$')

##
r_cmip5 = np.corrcoef(cmip6_deltaCstau_ssp585.flatten(), cmip6_deltaCsRh_ssp585.flatten())[0, 1]**2
print('CMIP6 delta', r_cmip5)
r_cmip5 = np.corrcoef(cmip6_fractionaldeltaCstau_ssp585.flatten(), cmip6_fractionaldeltaCsRh_ssp585.flatten())[0, 1]**2
print('CMIP6 frac', r_cmip5)    

        
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
leg1 = ax.legend(handles, labels, loc='center right', borderaxespad=0.2, bbox_to_anchor=(1.8, 0.5), title='CMIP6', fontsize=40)
plt.gca().add_artist(leg1)

handles2 = []
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=30, color='k', label='SSP126')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='v', markersize=30, color='k', label='SSP245')])
handles2.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=30, color='k', label='SSP585')])
labels2 = ['SSP126', 'SSP245', 'SSP585']
leg2 = ax.legend(handles2, labels2, loc='upper right', borderaxespad=0.2, fontsize=40)
plt.gca().add_artist(leg2)


#%%
fig_figure1.savefig('figures/CsNPPVsCstau_cmip5cmip6_v1', bbox_inches='tight')
plt.close()
