#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2023
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots timeseries for dtau for each CMIP6 ESMs using the C4MIP simulations 1pctCO2 and 1pctCO2-bgc,
then compares the ratios between the to simulations.
(A measure of how much dtau is due to CO2, i.e. increased NPP input.)
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

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF
from rmv_cmip_analysis import annual_average, area_average, global_total, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(42,12))
gs = gspec.GridSpec(1, 3, figure=fig_figure1, wspace=0.3)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':7,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':38,
    'xtick.labelsize':38,
    'ytick.labelsize':38,
    'font.size':38,
}
plt.rcParams.update(params)


#%%

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)

# C4MIP simulations
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc'] 


# CO2 array
co2 = np.zeros((150))
for i in range(0, 150):
        if i == 0:
             co2[i] = 285
        else:
            co2[i] = co2[i-1]*1.01

# global region used for global averages
region_global = [0, 360, -90,  90]

fp_percentages = np.zeros((n_models))


#%%
# for loop for each CMIP6 model
for model_i in range(0, n_models):
        model = cmip6_models[model_i] # seleting the models

        # Loop through each C4MIP simulation
        for c4mip_option in range(0, len(C4MIP_simulation)):
            c4mip = C4MIP_simulation[c4mip_option] # selecting the rcp scenario

            print(model, c4mip)    
    
            #%%
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)

            
            # cSoil
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cSoil_Emon_'+model+'_'+c4mip+'*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = annual_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            # cLitter
            if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cLitter_Lmon_'+model+'_'+c4mip+'*', model)
                cLitter_cube = open_netCDF(cLitter_cube)
                cLitter_cube = annual_average(cLitter_cube)
                cLitter_cube = global_total_percentage(cLitter_cube, landfrac=landfraction, latlon_cons=None)
                cLitter_data = cLitter_cube.data
                #
                Cs_data = cSoil_data + cLitter_data
            else:
                Cs_data = cSoil_data.copy()
            #rolling mean
            Cs_data_rm = np.zeros((len(Cs_data)-4))
            for i in range(2, len(Cs_data)-2):
                Cs_average = np.mean(Cs_data[i-2:i+2])
                Cs_data_rm[i-2] = Cs_average


            # Time dimension
            if model == 'ACCESS-ESM1-5':
                time_dimension = cSoil_cube.coord('year').points
                time_dimension = time_dimension - 100
            
            # tau 
            rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/rh_Lmon_'+model+'_'+c4mip+'*', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = annual_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube, landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data*86400.*360.
            #rolling mean
            rh_data_rm = np.zeros((len(rh_data)-4))
            for i in range(2, len(rh_data)-2):
                rh_average = np.mean(rh_data[i-2:i+2])
                rh_data_rm[i-2] = rh_average
            
            # dtau
            tau_data = Cs_data_rm/rh_data_rm
            delta_tau = tau_data - np.mean(tau_data[0:5])
            
            
            ax = fig_figure1.add_subplot(gs[0, c4mip_option])
            
            time_dimension = time_dimension[0:135]
            delta_tau = delta_tau[0:135]
            
            ax.plot(time_dimension, delta_tau, color=model_colors[model_i])
            
            
            plt.xlabel('Time (yr)')
            plt.xlim((0,135))
            plt.ylim((-32,5))
            
            
            if c4mip == '1pctCO2':
                dtau_1pct = delta_tau.copy()
                ax.set_title('(a) 1pctCO2', y=1.1, fontweight='bold')
                plt.ylabel(r'$\Delta \tau_{s}$ (yr)')
            else:
                dtau_bgc = delta_tau.copy()
                ax.set_title('(b) 1pctCO2-bgc', y=1.1, fontweight='bold')
                plt.ylabel(r'$\Delta \tau_{s}^{BGC}$ (yr)')
                
#            if model=='UKESM1-0-LL' and c4mip=='1pctCO2-bgc':
#                handles = []
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='peachpuff', label='ACCESS-ESM1-5')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#fb8072', label='BCC-CSM2-MR')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM5')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='dodgerblue', label='CESM2')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkcyan', label='GFDL-ESM4')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM6A-LR')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='olive', label='MIROC-ES2L')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='gold', label='MPI-ESM1-2-LR')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='orange', label='NorESM2-LM')])
#                handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label='UKESM1-0-LL')])
#                labels = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
#                leg1 = ax.legend(handles, labels, loc='lower center', ncol=2, borderaxespad=0.2, bbox_to_anchor=(-0.15, -0.75), title='CMIP6 ESMs', fontsize=38)
#                plt.gca().add_artist(leg1)

            
        #%%
        
        ax = fig_figure1.add_subplot(gs[0, 2])
        ratio = dtau_bgc[-1]/dtau_1pct[-1]
#        ratio = np.mean(dtau_bgc[130:135])/np.mean(dtau_1pct[130:135])
#        ax.plot(time_dimension, ratio, color=model_colors[model_i])
        print(model, ratio*100)
        
        
        fp_percentages[model_i] = ratio#*100#np.array([70, 79, 53, 89, 58, 60, 50, 63, 93, 75])
        
        
label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']

x = np.arange(10)  # the label locations
width = 0.5  # the width of the bars

rects = plt.bar(x, fp_percentages, width, color=['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen'])#color='cornflowerblue')


#ax.set_ylabel(r'$\Delta \tau_{s}^{BGC}$ / $\Delta \tau_{s}$')
ax.set_ylabel(r'$\frac{\Delta \tau_{s}^{BGC}}{\Delta \tau_{s}}$')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
plt.ylim((0, 1))
#plt.set_title('SSP126', fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')

ax.set_title(r'(c) Fraction of $\Delta \tau_{s}$ due to $\Delta CO_{2}$', y=1.1, fontweight='bold')
        
#%%
fig_figure1.savefig('figures/fp_contribution_c4mip_v5', bbox_inches='tight')
plt.close()