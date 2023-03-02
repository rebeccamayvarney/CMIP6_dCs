#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 2021
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dnpp timseries for CMIP5 and CMIP6 ESMs in future SSP and RCP simulations.
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
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_cmip5, combine_netCDF_model, open_netCDF, numpy_to_cube
from rmv_cmip_analysis import select_time, time_average, annual_average, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(48,28))
gs = gspec.GridSpec(2, 3, figure=fig_figure1, hspace=0.4, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':7,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':58,
    'xtick.labelsize':58,
    'ytick.labelsize':58,
    'font.size':58,
}
plt.rcParams.update(params)


# global region used for global averages
region_global = [0, 360, -90,  90]


#%% CMIP5

# rcp scenarios
rcp_options = ['rcp26', 'rcp45', 'rcp85']
rcp_options_length = len(rcp_options)
rcp_titles = ['RCP2.6', 'RCP4.5', 'RCP8.5']

# models
cmip5_models = cmip5_models = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)
model_colors_cmip5 = ['darkblue', '#80b1d3', 'darkcyan', '#8dd3c7', 'darkseagreen', 'darkgreen', 'olive', 'gold', 'orange']
        
                      
row = 0
# Loop through each rcp run being considered
for rcp_option in range(0, rcp_options_length):
    rcp = rcp_options[rcp_option] # selecting the rcp scenario

    # subplot    
    ax = fig_figure1.add_subplot(gs[row, rcp_option])

    # for loop for each CMIP5 model
    for model_i in range(0, n_models):
            model = cmip5_models[model_i] # seleting the models
    
            print(rcp, model)
            if rcp=='rcp26':
                if model=='MPI-ESM-LR':
                    continue
                
            # NPP (Historical) 
            npp_historical_cube_new = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/npp_Lmon_'+model+'_historical*', 'net_primary_productivity_of_carbon', model)
            npp_historical_cube_new = open_netCDF(npp_historical_cube_new)
            npp_historical_cube_new = select_time(npp_historical_cube_new, 1995, 2005)
            npp_historical_time_av_cube_new = time_average(npp_historical_cube_new)
            npp_historical_time_av_data_new = npp_historical_time_av_cube_new.data*86400.*360.
            # NPP (Future)
            npp_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/npp_Lmon_'+model+'_'+rcp+'_*', 'net_primary_productivity_of_carbon', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = annual_average(npp_cube)
            npp_time_av_data = npp_cube.data*86400.*360.
            
            ##
            time_dimension = npp_cube.coord('year').points
            cube_save = npp_cube.copy()
            
    
            #%% Finding the timeseries data
    
            # Modelled delta npp
            delta_npp_actual_cmip5 = npp_time_av_data - npp_historical_time_av_data_new
            # Masking invalid values
            delta_npp_actual_cmip5 = np.ma.masked_invalid(delta_npp_actual_cmip5)
            # convert numpy array to cube
            delta_npp_actual_cmip5_cube = numpy_to_cube(delta_npp_actual_cmip5, cube_save, 3)
            # Calculating the global averaged values delta Cs
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip5_data/sftlf_fx_'+model+'_*', model)
            actual_delta_npp_global_cmip5_cube = global_total_percentage(delta_npp_actual_cmip5_cube, landfrac=landfraction, latlon_cons=None)
            actual_delta_npp_global_cmip5_data = actual_delta_npp_global_cmip5_cube.data
            
            
            #%% Plotting
            
            if model=='CCSM4' or model=='NorESM1-M':
                plt.plot(time_dimension, actual_delta_npp_global_cmip5_data, color=model_colors_cmip5[model_i], linestyle='dashed')
            else:
                plt.plot(time_dimension, actual_delta_npp_global_cmip5_data, color=model_colors_cmip5[model_i])
            
            plt.title(rcp_titles[rcp_option], fontweight='bold')
            plt.xlabel('Time (yr)')
            plt.ylabel(r'$\Delta$ NPP (PgC yr$^{-1}$)')
            plt.xlim((2015,2100))
            plt.ylim((-5,75))
        
##
handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkblue', label='BNU-ESM')])
#handles.extend([Line2D([0,0],[0,0], linewidth=20, color='dodgerblue', label='CCSM4')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM2')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkcyan', label='GFDL-ESM2G')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#8dd3c7', label='GISS-E2-R')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label='HadGEM2-ES')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM5A-LR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='olive', label='MIROC-ESM')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='gold', label='MPI-ESM-LR')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='orange', label='NorESM1-M')])
labels = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
leg1 = ax.legend(handles, labels, handlelength=3, loc='center right', borderaxespad=0.2, bbox_to_anchor=(2.15, 0.55), title='CMIP5', fontsize=58)
plt.gca().add_artist(leg1)


#%% CMIP6

# ssp scenarios
ssp_options = ['ssp126', 'ssp245', 'ssp585']
ssp_options_length = len(ssp_options)
ssp_titles = ['SSP126', 'SSP245', 'SSP585']

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'red', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


row = 1
# Loop through each ssp run being considered
for ssp_option in range(0, ssp_options_length):
    ssp = ssp_options[ssp_option] # selecting the ssp scenario

    # subplot    
    ax = fig_figure1.add_subplot(gs[row, ssp_option])

    # for loop for each CMIP5 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i] # seleting the models
    
            print(ssp, model)
            if ssp=='ssp245' or ssp=='ssp585':
                if model=='GFDL-ESM4':
                    continue
                
            # NPP (Historical)
            npp_historical_cube_new = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
            npp_historical_cube_new = open_netCDF(npp_historical_cube_new)
            npp_historical_cube_new = select_time(npp_historical_cube_new, 2000, 2015)
            npp_historical_time_av_cube_new = time_average(npp_historical_cube_new)
            npp_historical_time_av_data_new = npp_historical_time_av_cube_new.data*86400.*360.
            # NPP (Future)
            npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_'+ssp+'_*', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = annual_average(npp_cube)
            npp_time_av_data = npp_cube.data*86400.*360.
            
            ##
            time_dimension = npp_cube.coord('year').points
            cube_save = npp_cube.copy()
            
    
            #%% Finding the timeseries data
    
            # Modelled delta npp
            delta_npp_actual_cmip6 = npp_time_av_data - npp_historical_time_av_data_new
            # Masking invalid values
            delta_npp_actual_cmip6 = np.ma.masked_invalid(delta_npp_actual_cmip6)
            # convert numpy array to cube
            delta_npp_actual_cmip6_cube = numpy_to_cube(delta_npp_actual_cmip6, cube_save, 3)
            # Calculating the global averaged values delta Cs
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)
            actual_delta_npp_global_cmip6_cube = global_total_percentage(delta_npp_actual_cmip6_cube, landfrac=landfraction, latlon_cons=None)
            actual_delta_npp_global_cmip6_data = actual_delta_npp_global_cmip6_cube.data
            
            
            
#            # PERCENTAGE INCREASE
#            npp_historical_av = np.ma.masked_invalid(npp_historical_time_av_data_new)
#            npp_historical_av_cube = numpy_to_cube(npp_historical_av, cube_save, 2)
#            npp_historical_av_cube = global_total_percentage(npp_historical_av_cube, landfrac=landfraction, latlon_cons=None)
#            npp_historical_av_data = npp_historical_av_cube.data
#            #
#            npp_future_av = np.ma.masked_invalid(npp_time_av_data)
#            npp_future_av_cube = numpy_to_cube(npp_future_av, cube_save, 3)
#            npp_future_av_cube = global_total_percentage(npp_future_av_cube, landfrac=landfraction, latlon_cons=None)
#            npp_future_av_data2 = npp_future_av_cube.data
#            #
#            npp_percent_change = (npp_future_av_data2[-1] - npp_historical_av_data[0]) / npp_future_av_data2[-1]
#            print(model, npp_percent_change*100)
            
            
            #%% Plotting
            
            if model=='ACCESS-ESM1-5' or model=='CESM2' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='UKESM1-0-LL':
                plt.plot(time_dimension, actual_delta_npp_global_cmip6_data, color=model_colors[model_i], linestyle='dashed')
            else:
                plt.plot(time_dimension, actual_delta_npp_global_cmip6_data, color=model_colors[model_i])
                
            plt.title(ssp_titles[ssp_option], fontweight='bold')
            plt.xlabel('Time (yr)')
            plt.ylabel(r'$\Delta$ NPP (PgC yr$^{-1}$)')
            plt.xlim((2015,2100))
            plt.ylim((-5,75))
            
            if ssp=='ssp245':
                handles2 = []
                handles2.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=5, color='k', label='Nitrogen cycle')])
                labels2 = ['Nitrogen cycle']
                leg2 = ax.legend(handles2, labels2, loc='lower center', borderaxespad=0.2, bbox_to_anchor=(0.5, -0.5), fontsize=58)
                plt.gca().add_artist(leg2)
        
        
## 
handles = []
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='peachpuff', label='ACCESS-ESM1-5')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#fb8072', label='BCC-CSM2-MR')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM5')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='dodgerblue', label='CESM2')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='red', label='CNRM-ESM2-1')])
#handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkcyan', label='GFDL-ESM4')])
handles.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM6A-LR')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='olive', label='MIROC-ES2L')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='gold', label='MPI-ESM1-2-LR')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='orange', label='NorESM2-LM')])
handles.extend([Line2D([0,0],[0,0], linestyle='dashed', linewidth=20, color='darkseagreen', label='UKESM1-0-LL')])
labels = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
leg1 = ax.legend(handles, labels, handlelength=3, loc='center right', borderaxespad=0.2, bbox_to_anchor=(2.2, 0.45), title='CMIP6', fontsize=58)
plt.gca().add_artist(leg1)

        
#%%
fig_figure1.savefig('figures/dNPP_timeseries_cmip5cmip6_v1', bbox_inches='tight')
plt.close()