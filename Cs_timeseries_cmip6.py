#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on April 2021
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots mean dCs in future SSP runs for CMIP6 and CMIP5 (RCP) ESMs.
"

"""

#%%

# Analysis imports
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF
from rmv_cmip_analysis import annual_average, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%%

# Set up subplot figure
fig = plt.figure(1, figsize=(28,48))
fig.tight_layout()
gs = gspec.GridSpec(5, 2, figure=fig, hspace=0.6, wspace=0.5)
column_1 = 0
row_1 = 0
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':10,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':52,
    'xtick.labelsize':52,
    'ytick.labelsize':52,
    'font.size':52,
#    'text.usetex': False,
#    "svg.fonttype": 'none'
}
plt.rcParams.update(params)


#%% CMIP6

# models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

# ssp scenarios
ssp_options = ['ssp126', 'ssp245', 'ssp585']
ssp_options_length = len(ssp_options)


# for loop for each CMIP^ model
for model_i in range(0, n_models):
    model = cmip6_models[model_i] # seleting the models
    #
    ax = fig.add_subplot(gs[row_1, column_1])

    # Loop through each ssp run being considered
    for ssp_option in range(0, ssp_options_length):
        ssp = ssp_options[ssp_option] # selecting the ssp scenario

        print(ssp, model)
        if ssp=='ssp245' or ssp=='ssp585':
            if model=='GFDL-ESM4':
                continue
        
        # model land fraction
        landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)

        # Soil Carbon (cSoil)
        cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_'+ssp+'_*', model)
        cSoil_cube = open_netCDF(cSoil_cube)
        cSoil_cube = annual_average(cSoil_cube)
        time_dimension = cSoil_cube.coord('year').points
        cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
        cSoil_data = cSoil_cube.data

        # Litter Carbon (cLitter)
        if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM':
            cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_'+ssp+'*', model)
            cLitter_cube = open_netCDF(cLitter_cube)        
            cLitter_cube = annual_average(cLitter_cube)
            cLitter_cube = global_total_percentage(cLitter_cube, landfrac=landfraction, latlon_cons=None)
            cLitter_data = cLitter_cube.data
            #
            soil_carbon_future_data = cSoil_data + cLitter_data
        else:
            soil_carbon_future_data = cSoil_data.copy()


        #%% Plotting
        if ssp == 'ssp126':
            ax.plot(time_dimension, soil_carbon_future_data, color='darkkhaki')
        if ssp == 'ssp245':
            ax.plot(time_dimension, soil_carbon_future_data, color='mediumseagreen')
        if ssp == 'ssp585':
            ax.plot(time_dimension, soil_carbon_future_data, color='darkgreen')
        ax.set_title(model)
        
    ax.set_ylabel(r'$C_{\mathrm{s}}$ (PgC)')
    ax.set_xlabel('Time (yr)')
    #plt.ylim((0, ))
    plt.xlim((2015,2100))

    #
    row_1 += 1 
    if row_1==5: 
        column_1 += 1
        row_1 = 0

#%%
        
handle = []
handle.extend([Line2D([0,0],[0,0], linewidth=20, color='darkkhaki', label='SSP126')])
handle.extend([Line2D([0,0],[0,0], linewidth=20, color='mediumseagreen', label='SSP245')])
handle.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='SSP585')])
label = ['SSP126', 'SSP245', 'SSP585']
leg = plt.legend(handle, label, loc='lower center', ncol=3, bbox_to_anchor=(-0.16, -0.95), fontsize=52)
plt.gca().add_artist(leg)

fig.savefig('figures/Cs_timeseries_cmip5cmip6_v1', bbox_inches='tight')
plt.close()
