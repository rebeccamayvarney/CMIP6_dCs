#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dCs,npp Vs dCs,tau for each CMIP6 ESMs using the C4MIP simulations 1pctCO2, 1pctCO2-bgc, 1pctCO2-rad.
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


#%% Loading the data

# 1pctCO2 (2xCO2)
cmip6_deltaCs_1pctCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_2xCO2.npy')
cmip6_deltaCstau_1pctCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_2xCO2.npy')
# 1pctCO2-bgc (2xCO2)
cmip6_deltaCs_1pctCO2_bgc = np.load('saved_data/cmip6_Cs_1pctCO2_bgc_2xCO2.npy')
cmip6_deltaCstau_1pctCO2_bgc = np.load('saved_data/cmip6_Cstau_1pctCO2_bgc_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_bgc = np.load('saved_data/cmip6_CsNPP_1pctCO2_bgc_2xCO2.npy')
# 1pctCO2-rad (2xCO2)
cmip6_deltaCs_1pctCO2_rad = np.load('saved_data/cmip6_Cs_1pctCO2_rad_2xCO2.npy')
cmip6_deltaCstau_1pctCO2_rad = np.load('saved_data/cmip6_Cstau_1pctCO2_rad_2xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_rad = np.load('saved_data/cmip6_CsNPP_1pctCO2_rad_2xCO2.npy')

# 1pctCO2 (4xCO2)
cmip6_deltaCs_1pctCO2_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_4xCO2.npy')
# 1pctCO2-bgc (4xCO2)
cmip6_deltaCs_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_bgc_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_bgc_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_bgc_4xCO2.npy')
# 1pctCO2-rad (4xCO2)
cmip6_deltaCs_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_Cs_1pctCO2_rad_4xCO2.npy')
cmip6_deltaCstau_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_Cstau_1pctCO2_rad_4xCO2.npy')
cmip6_deltaCsNPP_1pctCO2_rad_4xCO2 = np.load('saved_data/cmip6_CsNPP_1pctCO2_rad_4xCO2.npy')


#%%
# Setting up the figure

fig_figure1 = plt.figure(1, figsize=(50,56))
gs = gspec.GridSpec(3, 2, figure=fig_figure1, hspace=0.4, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':84,
    'xtick.labelsize':84,
    'ytick.labelsize':84,
    'font.size':84,
}
plt.rcParams.update(params)


#%%

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors_cmip6 = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


# for loop for each CMIP6 model
for model_i in range(0, n_models):
    model = cmip6_models[model_i]


    #%% 1pctCO2 (2xCO2)
    ax = fig_figure1.add_subplot(gs[0, 0])
    
    plt.plot(cmip6_deltaCstau_1pctCO2[model_i], cmip6_deltaCsNPP_1pctCO2[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}$ (PgC)')
        plt.title(r'(a) 2xCO$_{2}$', y=1.15, fontweight='bold')
        
        # linear regression
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2, cmip6_deltaCsNPP_1pctCO2, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)
    
    
    #%% 1pctCO2-bgc (2xCO2)
    ax = fig_figure1.add_subplot(gs[1, 0])
    
    plt.plot(cmip6_deltaCstau_1pctCO2_bgc[model_i], cmip6_deltaCsNPP_1pctCO2_bgc[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}^{BGC}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}^{BGC}$ (PgC)')
        
        # linear regressions
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2_bgc, cmip6_deltaCsNPP_1pctCO2_bgc, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2_bgc)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)
    
    
    #%% 1pctCO2-rad (2xCO2)
    ax = fig_figure1.add_subplot(gs[2, 0])
    
    plt.plot(cmip6_deltaCstau_1pctCO2_rad[model_i], cmip6_deltaCsNPP_1pctCO2_rad[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}^{RAD}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}^{RAD}$ (PgC)')
        
        # linear reg
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2_rad, cmip6_deltaCsNPP_1pctCO2_rad, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2_rad)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)
        

    #%% 1pctCO2 (4xCO2)
    ax = fig_figure1.add_subplot(gs[0, 1])
    
    plt.plot(cmip6_deltaCstau_1pctCO2_4xCO2[model_i], cmip6_deltaCsNPP_1pctCO2_4xCO2[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}$ (PgC)')
        plt.title(r'(b) 4xCO$_{2}$', y=1.15, fontweight='bold')
        
        # linear regression
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2_4xCO2, cmip6_deltaCsNPP_1pctCO2_4xCO2, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2_4xCO2)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)
    
    
    #%% 1pctCO2-bgc (4xCO2)
    ax = fig_figure1.add_subplot(gs[1, 1])
    
    plt.plot(cmip6_deltaCstau_1pctCO2_bgc_4xCO2[model_i], cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}^{BGC}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}^{BGC}$ (PgC)')
        
        # linear regression
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2_bgc_4xCO2, cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2_bgc_4xCO2)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)
    
    
    #%% 1pctCO2-rad (4xCO2)
    ax = fig_figure1.add_subplot(gs[2, 1])
    
    plt.plot(cmip6_deltaCstau_1pctCO2_rad_4xCO2[model_i], cmip6_deltaCsNPP_1pctCO2_rad_4xCO2[model_i], color=model_colors_cmip6[model_i], marker='o', markersize=40)
    
    if model=='UKESM1-0-LL':   
        plt.xlabel(r'$\Delta C_{s, \tau}^{RAD}$ (PgC)')
        plt.ylabel(r'$\Delta C_{s, NPP}^{RAD}$ (PgC)')
        
        # linear regression
        m_1 = np.polyfit(cmip6_deltaCstau_1pctCO2_rad_4xCO2, cmip6_deltaCsNPP_1pctCO2_rad_4xCO2, 1)
        m = m_1[0]
        print(m)
        m_function = np.poly1d(m_1)
        sorted_tas = np.sort(cmip6_deltaCstau_1pctCO2_rad_4xCO2)
        #
        plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=5)


#%%
r_cmip6 = np.corrcoef(cmip6_deltaCstau_1pctCO2_4xCO2, cmip6_deltaCsNPP_1pctCO2_4xCO2)[0, 1]**2
print('full', r_cmip6)
r_cmip6 = np.corrcoef(cmip6_deltaCstau_1pctCO2_bgc_4xCO2, cmip6_deltaCsNPP_1pctCO2_bgc_4xCO2)[0, 1]**2
print('BGC', r_cmip6)
r_cmip6 = np.corrcoef(cmip6_deltaCstau_1pctCO2_rad_4xCO2, cmip6_deltaCsNPP_1pctCO2_rad_4xCO2)[0, 1]**2
print('RAD', r_cmip6)
        

#%%

handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='peachpuff', label='ACCESS-ESM1-5')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='#fb8072', label='BCC-CSM2-MR')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='#80b1d3', label='CanESM5')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='dodgerblue', label='CESM2')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='darkcyan', label='GFDL-ESM4')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='darkgreen', label='IPSL-CM6A-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='olive', label='MIROC-ES2L')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='gold', label='MPI-ESM1-2-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='orange', label='NorESM2-LM')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='darkseagreen', label='UKESM1-0-LL')])
labels = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=5, borderaxespad=0.2, bbox_to_anchor=(-0.25, -0.95), title='CMIP6', fontsize=84)
plt.gca().add_artist(leg1)

        
fig_figure1.savefig('figures/CsNPP_Vs_Cstau_C4MIP_v2', bbox_inches='tight')
plt.close()