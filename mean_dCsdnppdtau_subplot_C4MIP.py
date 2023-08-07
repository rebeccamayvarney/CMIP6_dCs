#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots timeseries for dCs, dNPP, dtau for each CMIP6 ESMs using the C4MIP simulations 1pctCO2, 1pctCO2-bgc, 1pctCO2-rad.
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

fig_figure1 = plt.figure(1, figsize=(48,32))
gs = gspec.GridSpec(3, 3, figure=fig_figure1, hspace=0.35, wspace=0.35)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':7,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':56,
    'xtick.labelsize':56,
    'ytick.labelsize':56,
    'font.size':56,
}
plt.rcParams.update(params)


#%%

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

# C4MIP simulations
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad'] 


# CO2 array
co2 = np.zeros((150))
for i in range(0, 150):
        if i == 0:
             co2[i] = 285
        else:
            co2[i] = co2[i-1]*1.01

# global region used for global averages
region_global = [0, 360, -90,  90]


#%%
# Loop through each C4MIP simulation
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option] # selecting the rcp scenario
    
    dCs_array_1 = np.zeros((n_models, 150))
    dnpp_array_1 = np.zeros((n_models, 150))
    dtau_array_1 = np.zeros((n_models, 150))
    
    dCs_array_2 = np.zeros((n_models, 150))
    dnpp_array_2 = np.zeros((n_models, 150))
    dtau_array_2 = np.zeros((n_models, 150))
    
    dCs_array_3 = np.zeros((n_models, 150))
    dnpp_array_3 = np.zeros((n_models, 150))
    dtau_array_3 = np.zeros((n_models, 150))


    # for loop for each CMIP6 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i] # seleting the models
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
            # dCs
            deltaCs = Cs_data - Cs_data[0]


            # Time dimension
            if model == 'ACCESS-ESM1-5':
                time_dimension = cSoil_cube.coord('year').points
                time_dimension = time_dimension - 100

            
            # npp
            npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/npp_Lmon_'+model+'_'+c4mip+'*', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = annual_average(npp_cube)
            npp_cube = global_total_percentage(npp_cube, landfrac=landfraction, latlon_cons=None)
            npp_data = npp_cube.data*86400.*360.
            # dnpp
            delta_npp = npp_data - npp_data[0]
            
            
            # tau 
            rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/rh_Lmon_'+model+'_'+c4mip+'*', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = annual_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube, landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data*86400.*360.
            # dtau
            tau_data = Cs_data/rh_data
            delta_tau = tau_data - tau_data[0]
            
            
            # Temperature
            tas_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/tas_Amon_'+model+'_'+c4mip+'*', model)
            tas_cube = open_netCDF(tas_cube)
            tas_cube = annual_average(tas_cube)
            tas_cube = area_average(tas_cube, [0, 360, -90,  90])
            tas_data = tas_cube.data - 273.15
            #dT
            deltaT = tas_data - tas_data[0]
            
            
            #%%
            
            if model=='CanESM5':
                deltaCs = deltaCs[:-1]
                delta_npp = delta_npp[:-1]
                delta_tau = delta_tau[:-1]
                
            if model=='BCC-CSM2-MR' or model=='CanESM5':
                deltaT = deltaT[:-1]
            
            if model=='MPI-ESM1-2-LR' and c4mip=='1pctCO2':
                deltaCs = deltaCs[:-15]
                delta_npp = delta_npp[:-15]
                delta_tau = delta_tau[:-15]
                
#            if model=='MPI-ESM1-2-LR' and c4mip=='1pctCO2-bgc':
#                time_dimension1 = time_dimension[:-10]
#                co2_1 = co2[:-10]
#                ax = fig_figure1.add_subplot(gs[0, c4mip_option])
#                plt.plot(co2_1, deltaCs, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[1, c4mip_option])
#                plt.plot(co2_1, delta_npp, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[2, c4mip_option])
#                plt.plot(co2_1, delta_tau, color=model_colors[model_i])#, linestyle='dashed')
#                continue
                
            if model=='MPI-ESM1-2-LR' and c4mip=='1pctCO2-rad':
                for i in range(10):
                    deltaCs.append(deltaCs[-1])
                dCs_array_3[model_i, :] = deltaCs
                dnpp_array_3[model_i, :] = delta_npp
                dtau_array_3[model_i, :] = delta_tau
#                x_axis_new = deltaT
#                x_axis_new_av = np.mean(x_axis_new.reshape(-1, 5), axis=1)
#                x_axis_new = np.repeat(x_axis_new_av, 5)
#                #
#                deltaCs_avg = np.mean(deltaCs.reshape(-1, 5), axis=1)
#                deltaCs = np.repeat(deltaCs_avg, 5)
#                delta_npp_avg = np.mean(delta_npp.reshape(-1, 5), axis=1)
#                delta_npp = np.repeat(delta_npp_avg, 5)
#                delta_tau_avg = np.mean(delta_tau.reshape(-1, 5), axis=1)
#                delta_tau = np.repeat(delta_tau_avg, 5)
#                
#                time_dimension1 = time_dimension[:-10]
#                ax = fig_figure1.add_subplot(gs[0, c4mip_option])
#                plt.plot(x_axis_new, deltaCs, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[1, c4mip_option])
#                plt.plot(x_axis_new, delta_npp, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[2, c4mip_option])
#                plt.plot(x_axis_new, delta_tau, color=model_colors[model_i])#, linestyle='dashed')
#                continue
                
#            if model=='NorESM2-LM' and c4mip=='1pctCO2-bgc':
#                time_dimension1 = time_dimension[:-10]
#                co2_1 = co2[:-10]
#                ax = fig_figure1.add_subplot(gs[0, c4mip_option])
#                plt.plot(co2_1, deltaCs, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[1, c4mip_option])
#                plt.plot(co2_1, delta_npp, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[2, c4mip_option])
#                plt.plot(co2_1, delta_tau, color=model_colors[model_i])#, linestyle='dashed')
#                continue
#                
            if model=='NorESM2-LM' and c4mip=='1pctCO2-rad':
#                x_axis_new = deltaT
#                x_axis_new_av = np.mean(x_axis_new.reshape(-1, 10), axis=1)
#                x_axis_new = np.repeat(x_axis_new_av, 10)
#                #
#                deltaCs_avg = np.mean(deltaCs.reshape(-1, 10), axis=1)
#                deltaCs = np.repeat(deltaCs_avg, 10)
#                delta_npp_avg = np.mean(delta_npp.reshape(-1, 10), axis=1)
#                delta_npp = np.repeat(delta_npp_avg, 10)
#                delta_tau_avg = np.mean(delta_tau.reshape(-1, 10), axis=1)
#                delta_tau = np.repeat(delta_tau_avg, 10)
#                
#                ax = fig_figure1.add_subplot(gs[0, c4mip_option])
#                plt.plot(x_axis_new, deltaCs, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[1, c4mip_option])
#                plt.plot(x_axis_new, delta_npp, color=model_colors[model_i])#, linestyle='dashed')
#                ax = fig_figure1.add_subplot(gs[2, c4mip_option])
#                plt.plot(x_axis_new, delta_tau, color=model_colors[model_i])#, linestyle='dashed')
#                continue


            #%% fill arrays
            
            if c4mip=='1pctCO2':
                dCs_array_1[model_i, :] = deltaCs
                dnpp_array_1[model_i, :] = delta_npp
                dtau_array_1[model_i, :] = delta_tau
            elif c4mip=='1pctCO2-bgc':
                dCs_array_2[model_i, :] = deltaCs
                dnpp_array_2[model_i, :] = delta_npp
                dtau_array_2[model_i, :] = delta_tau
            elif c4mip=='1pctCO2-rad':
                dCs_array_3[model_i, :] = deltaCs
                dnpp_array_3[model_i, :] = delta_npp
                dtau_array_3[model_i, :] = delta_tau
    
    if c4mip=='1pctCO2':
        mean_dCs1 = np.mean(dCs_array_1, axis=0)
        mean_dnpp1 = np.mean(dnpp_array_1, axis=0)
        mean_dtau1 = np.mean(dtau_array_1, axis=0)
        
        std_dCs1 = np.std(dCs_array_1, axis=0)
        std_dnpp1 = np.std(dnpp_array_1, axis=0)
        std_dtau1 = np.std(dtau_array_1, axis=0)
        
    elif c4mip=='1pctCO2-bgc':
        mean_dCs2 = np.mean(dCs_array_2, axis=0)
        mean_dnpp2 = np.mean(dnpp_array_2, axis=0)
        mean_dtau2 = np.mean(dtau_array_2, axis=0)
        
        std_dCs2 = np.std(dCs_array_2, axis=0)
        std_dnpp2 = np.std(dnpp_array_2, axis=0)
        std_dtau2 = np.std(dtau_array_2, axis=0)
        
    elif c4mip=='1pctCO2-rad':
        mean_dCs3 = np.mean(dCs_array_3, axis=0)
        mean_dnpp3 = np.mean(dnpp_array_3, axis=0)
        mean_dtau3 = np.mean(dtau_array_3, axis=0)
        
        std_dCs3 = np.std(dCs_array_3, axis=0)
        std_dnpp3 = np.std(dnpp_array_3, axis=0)
        std_dtau3 = np.std(dtau_array_3, axis=0)
    
    
    #%% dCs
    ax = fig_figure1.add_subplot(gs[0, c4mip_option])
    
    if c4mip=='1pctCO2':
        x_axis = time_dimension
    elif c4mip=='1pctCO2-bgc':
        x_axis = co2
    elif c4mip=='1pctCO2-rad':
        x_axis = deltaT
#                x_axis_av = np.mean(x_axis.reshape(-1, 5), axis=1)
#                x_axis = np.repeat(x_axis_av, 5)
#                #
#                deltaCs_avg = np.mean(deltaCs.reshape(-1, 5), axis=1)
#                deltaCs = np.repeat(deltaCs_avg, 5)
#                delta_npp_avg = np.mean(delta_npp.reshape(-1, 5), axis=1)
#                delta_npp = np.repeat(delta_npp_avg, 5)
#                delta_tau_avg = np.mean(delta_tau.reshape(-1, 5), axis=1)
#                delta_tau = np.repeat(delta_tau_avg, 5)

    
    if c4mip=='1pctCO2':
        plt.plot(x_axis, mean_dCs1, color='b')
        plt.fill_between(mean_dCs1, mean_dCs1-std_dCs1, mean_dCs1+std_dCs1, color='b', alpha=0.25)
        ax.set_title('(a) 1pctCO2', y=1.15, fontweight='bold')
        plt.ylim((-25,700))
    if c4mip=='1pctCO2-bgc':
        plt.plot(x_axis, mean_dCs2, color='b')
        plt.fill_between(mean_dCs2, mean_dCs2-std_dCs2, mean_dCs2+std_dCs2, color='b', alpha=0.25)
        ax.set_title('(b) 1pctCO2-bgc', y=1.15, fontweight='bold')
        plt.ylim((-25,700))
    if c4mip=='1pctCO2-rad':
        plt.plot(x_axis, mean_dCs3, color='b')
        plt.fill_between(mean_dCs3, mean_dCs3-std_dCs3, mean_dCs3+std_dCs3, color='b', alpha=0.25)
        ax.set_title('(c) 1pctCO2-rad', y=1.15, fontweight='bold')
        plt.ylim((-425,25))

    
    plt.ylabel(r'$\Delta C_{s}$ (PgC)')
    if c4mip=='1pctCO2':
        plt.xlabel('Time (yr)')
        plt.xlim((0,150))
    elif c4mip=='1pctCO2-bgc':
        plt.xlabel(r'CO$_{2}$ (ppm)')
        plt.xlim((285,1100))
    elif c4mip=='1pctCO2-rad':
        plt.xlabel(r'$\Delta$T ($^{\circ}$C)')
        plt.xlim((0,4))


    #%% dnpp
    ax = fig_figure1.add_subplot(gs[1, c4mip_option])

    
    if c4mip=='1pctCO2':
        #plt.plot(x_axis, mean_dnpp1, color='g')
        plt.fill_between(mean_dnpp1, mean_dnpp1-std_dnpp1, mean_dnpp1+std_dnpp1, color='g', alpha=0.25)
        plt.ylim((-5,90))
    if c4mip=='1pctCO2-bgc':
        #plt.plot(x_axis, mean_dnpp2, color='g')
        plt.fill_between(mean_dnpp2, mean_dnpp2-std_dnpp2, mean_dnpp2+std_dnpp2, color='g', alpha=0.25)
        plt.ylim((-5,90))
    if c4mip=='1pctCO2-rad':
        plt.plot(x_axis, mean_dnpp3, color='g')
        plt.fill_between(mean_dnpp3, mean_dnpp3-std_dnpp3, mean_dnpp3+std_dnpp3, color='g', alpha=0.25)
        plt.ylim((-22,5))
    
    plt.ylabel(r'$\Delta$NPP (PgC yr$^{-1}$)')
    if c4mip=='1pctCO2':
        plt.xlabel('Time (yr)')
        plt.xlim((0,150))
    elif c4mip=='1pctCO2-bgc':
        plt.xlabel(r'CO$_{2}$ (ppm)')
        plt.xlim((285,1100))
    elif c4mip=='1pctCO2-rad':
        plt.xlabel(r'$\Delta$T ($^{\circ}$C)')
        plt.xlim((0,4))
    
    
#    if c4mip=='1pctCO2-bgc':
#        handles = []
#        handles.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label='Mean')])
#        handles.extend([Line2D([0,0],[0,0], linewidth=20, color='#fb8072', label='Standard Deviation')])
#        labels = ['Ensemble mean', 'Ensemble standard deviation']
#        leg1 = ax.legend(handles, labels, loc='lower center', ncol=2, borderaxespad=0.2, bbox_to_anchor=(0.5, -2.4), fontsize=56)
#        plt.gca().add_artist(leg1)
    
    
    #%% dtau
    ax = fig_figure1.add_subplot(gs[2, c4mip_option])
    
    
    if c4mip=='1pctCO2':
        plt.plot(x_axis, mean_dtau1, color='r')
        plt.fill_between(mean_dtau1, mean_dtau1-std_dtau1, mean_dtau1+std_dtau1, color='r', alpha=0.25)
        plt.ylim((-32,5))
    if c4mip=='1pctCO2-bgc':
        plt.plot(x_axis, mean_dtau2, color='r')
        plt.fill_between(mean_dtau2, mean_dtau2-std_dtau2, mean_dtau2+std_dtau2, color='r', alpha=0.25)
        plt.ylim((-32,5))
    if c4mip=='1pctCO2-rad':
        plt.plot(x_axis, mean_dtau3, color='r')
        plt.fill_between(mean_dtau3, mean_dtau3-std_dtau3, mean_dtau3+std_dtau3, color='r', alpha=0.25)
        plt.ylim((-11,5))


    plt.ylabel(r'$\Delta \tau_{s}$ (yr)')
    if c4mip=='1pctCO2':
        plt.xlabel('Time (yr)')
        plt.xlim((0,150))
    elif c4mip=='1pctCO2-bgc':
        plt.xlabel(r'CO$_{2}$ (ppm)')
        plt.xlim((285,1100))
    elif c4mip=='1pctCO2-rad':
        plt.xlabel(r'$\Delta$T ($^{\circ}$C)')
        plt.xlim((0,4))



#%%
fig_figure1.savefig('figures/mean_CsNPPtau_timeseries_c4mip_v1', bbox_inches='tight')
plt.close()