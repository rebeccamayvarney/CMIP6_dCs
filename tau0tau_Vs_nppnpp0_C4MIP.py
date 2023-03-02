#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots tau0*tau Vs npp*npp0 for each CMIP6 ESMs using the C4MIP simulations 1pctCO2, 1pctCO2-bgc, 1pctCO2-rad.
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
from rmv_cmip_analysis import annual_average, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%%

# C4MIP simulations
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc' , '1pctCO2-rad']
#C4MIP_labels = ['1pctCO2', '1pctCO2_bgc', '1pctCO2_rad']

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


#%% Set up subplot figure

fig = plt.figure(1, figsize=(32,64))
gs = gspec.GridSpec(5, 2, figure=fig, width_ratios=[1, 1], hspace=0.7, wspace=0.6)
n = 10
column_1 = 0
row_1 = 0
n_columns_1 = 3
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True 
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':72,
    'xtick.labelsize':72,
    'ytick.labelsize':72,
    'font.size':72,
}
plt.rcParams.update(params)


# Loop through each C4MIP simulation being considered
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option]
    #C4MIP_label = C4MIP_labels[c4mip_option]
    
    column_1 = 0
    row_1 = 0
    n_columns_1 = 3
    

    # for loop for each CMIP6 model
    for model_i, a, in zip(range(n_models), range(n)):
            model = cmip6_models[model_i]
            print(model, c4mip)
            
            # land fraction
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
            #
            Cs_data_trimmed = Cs_data[0:140]
            
            
            # Rh
            rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/rh_Lmon_'+model+'_'+c4mip+'*', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = annual_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube*86400.*360., landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data
            rh_data_trimmed = rh_data[0:140]

            # tau            
            tau_data_trimmed = Cs_data_trimmed/rh_data_trimmed
            tau_data_0 = tau_data_trimmed[0]
            #
            tau_data_trimmed_1 = tau_data_0/tau_data_trimmed

            
            # NPP
            npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/npp_Lmon_'+model+'_'+c4mip+'*', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = annual_average(npp_cube)
            npp_cube = global_total_percentage(npp_cube*86400.*360., landfrac=landfraction, latlon_cons=None)
            npp_data = npp_cube.data
            npp_data_trimmed = npp_data[0:140]
            npp_data_0 = npp_data_trimmed[0]
            #
            npp_data_trimmed_1 = npp_data_trimmed/npp_data_0
            
            
            ##
            ax = fig.add_subplot(gs[row_1, column_1])
            
            if c4mip=='1pctCO2':
                ax.plot(tau_data_trimmed_1, npp_data_trimmed_1, 'b.', markersize=40)
            elif c4mip=='1pctCO2-bgc':
                ax.plot(tau_data_trimmed_1, npp_data_trimmed_1, 'g.', markersize=40)
            elif c4mip=='1pctCO2-rad':
                ax.plot(tau_data_trimmed_1, npp_data_trimmed_1, 'r.', markersize=40)
            
            
            r = np.corrcoef(tau_data_trimmed_1, npp_data_trimmed_1)[0, 1]**2
            print(c4mip, model, r)

            
            ax.set_ylabel(r'NPP / NPP$_{0}$')
            ax.set_xlabel(r'$\tau_{s, 0}$ / $\tau_{s}$')
            #plt.ylim((1,2.5))
            #plt.xlim((1,2.5))
            ax.set_title(model)

            
            #%% increase row and column 
            row_1 += 1 
            if row_1==5:
                column_1 += 1
                row_1 = 0
                
##
handles = []
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='b', label=r'1')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='g', label=r'2')])
handles.extend([Line2D([0,0],[0,0], linewidth=40, color='r', label=r'3')])
labels = [r'Full 1% CO$_{2}$', r'BGC 1% CO$_{2}$', r'RAD 1% CO$_{2}$']
leg1 = ax.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(-0.3, -1.1), fontsize=72)
plt.gca().add_artist(leg1)    
      
      
#%%
fig.savefig('figures/NPPVstau_c4mip_v2', bbox_inches='tight')
plt.close()