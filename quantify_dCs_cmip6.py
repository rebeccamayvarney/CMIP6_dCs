#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 2021
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script quantifies and saves dCs and the breakdown values (dCs,npp / dCs,tau / dnppdtau, dCs,nep / dCs,tau_nep / dnepdtau)
for CMIP6 ESMs for SSP simulations.
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
from rmv_cmip_analysis import select_time, time_average, global_total_percentage

            
#%% CMIP6

# ssp scenarios
ssp_options = ['ssp126', 'ssp245', 'ssp585']
ssp_options_length = len(ssp_options)

# models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)


#%%
# Loop through each ssp run being considered
for ssp_option in range(0, ssp_options_length):
    ssp = ssp_options[ssp_option] # selecting the ssp scenario

    # Empty arrays
    cmip6_deltaCs = np.zeros((len(cmip6_models)))
    cmip6_deltaCstau = np.zeros((len(cmip6_models)))
    cmip6_deltaCsRh = np.zeros((len(cmip6_models)))
    cmip6_deltadelta = np.zeros((len(cmip6_models)))
    cmip6_remaining = np.zeros((len(cmip6_models)))
    cmip6_NEPtau = np.zeros((len(cmip6_models)))
    cmip6_deltaNEPdeltatau = np.zeros((len(cmip6_models)))
    cmip6_NEPdeltatau = np.zeros((len(cmip6_models)))
    cmip6_percent = np.zeros((len(cmip6_models)))
    #
    cmip6_deltaCs_frac = np.zeros((len(cmip6_models)))
    cmip6_deltaCstau_frac = np.zeros((len(cmip6_models)))
    cmip6_deltaCsRh_frac = np.zeros((len(cmip6_models)))
    cmip6_deltadelta_frac = np.zeros((len(cmip6_models)))
    
    
    # for loop for each CMIP6 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i] # seleting the models
            print(ssp, model)

            # land fraction
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_*', model)
            
            
            #%%
            # cSoil (historical)
            cSoil_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
            cSoil_historical_cube = open_netCDF(cSoil_historical_cube)
            cSoil_historical_cube = select_time(cSoil_historical_cube, 1995, 2005)
            cSoil_historical_cube = time_average(cSoil_historical_cube)
            cSoil_historical_cube = global_total_percentage(cSoil_historical_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_historical_data = cSoil_historical_cube.data
            
            # cLitter (historical)
            if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_historical*', model)
                cLitter_historical_cube = open_netCDF(cLitter_historical_cube)     
                cLitter_historical_cube = select_time(cLitter_historical_cube, 1995, 2005)
                cLitter_historical_cube = time_average(cLitter_historical_cube)
                cLitter_historical_cube = global_total_percentage(cLitter_historical_cube, landfrac=landfraction, latlon_cons=None)
                cLitter_historical_data = cLitter_historical_cube.data
                #
                soil_carbon_historical_data = cSoil_historical_data + cLitter_historical_data
            else:
                soil_carbon_historical_data = cSoil_historical_data.copy()
            
            
            # Rh (historical)
            rh_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/rh_Lmon_'+model+'_historical*', model)
            rh_historical_cube = open_netCDF(rh_historical_cube)
            rh_historical_cube = select_time(rh_historical_cube, 1995, 2005)
            rh_historical_cube = time_average(rh_historical_cube)
            rh_historical_cube = global_total_percentage(rh_historical_cube, landfrac=landfraction, latlon_cons=None)
            rh_historical_data = rh_historical_cube.data*86400.*360.
            #
            tau_0_data = soil_carbon_historical_data/rh_historical_data
            
            # NPP (historical)
            npp_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/npp_Lmon_'+model+'_historical*', model)
            npp_historical_cube = open_netCDF(npp_historical_cube)
            npp_historical_cube = select_time(npp_historical_cube, 1995, 2005)
            npp_historical_cube = time_average(npp_historical_cube)
            npp_historical_cube = global_total_percentage(npp_historical_cube, landfrac=landfraction, latlon_cons=None)
            npp_historical_data = npp_historical_cube.data*86400.*360.
            
            
            # cSoil (future)
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_historical_'+ssp+'/cSoil_Emon_'+model+'_*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = select_time(cSoil_cube, 2090, 2100)
            cSoil_cube = time_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            
            # cLitter (future)
            if  model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_historical_'+ssp+'/cLitter_Lmon_'+model+'_*', model)
                cLitter_cube = open_netCDF(cLitter_cube)     
                cLitter_cube = select_time(cLitter_cube, 2090, 2100)
                cLitter_cube = time_average(cLitter_cube)
                cLitter_cube = global_total_percentage(cLitter_cube, landfrac=landfraction, latlon_cons=None)
                cLitter_data = cLitter_cube.data
                #
                soil_carbon_future_data = cSoil_data + cLitter_data
            else:
                soil_carbon_future_data = cSoil_data.copy()
                
            # Rh (future)
            rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_historical_'+ssp+'/rh_Lmon_'+model+'_*', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = select_time(rh_cube, 2090, 2100)
            rh_cube = time_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube, landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data*86400.*360.
            #
            tau_data = soil_carbon_future_data/rh_data
            
            # NPP (future)
            npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_historical_'+ssp+'/npp_Lmon_'+model+'_*', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = select_time(npp_cube, 2090, 2100)
            npp_cube = time_average(npp_cube)
            npp_cube = global_total_percentage(npp_cube, landfrac=landfraction, latlon_cons=None)
            npp_data = npp_cube.data*86400.*360.
            
            
            #%%

            # dCs
            delta_Cs = soil_carbon_future_data - soil_carbon_historical_data
            
            # dCs,npp
            delta_rh = npp_data - npp_historical_data
            deltaCs_rh = delta_rh*tau_0_data
            
            # dCs,tau
            delta_tau = tau_data - tau_0_data
            deltaCs_tau = delta_tau*npp_historical_data
            
            # dnpp*dtau
            deltadelta = delta_rh*delta_tau
            
            # dNEP
            NEP_t = npp_data - rh_data
            NEP_0 = npp_historical_data - rh_historical_data
            # dCs,nep
            deltaCs_NEPtau = (NEP_t - NEP_0)*tau_0_data
            # dCs,tau_nep
            deltaCs_NEPdeltatau = NEP_0*delta_tau
            # dnep*dtau
            deltaNEPdeltatau = (NEP_t - NEP_0)*delta_tau
        
                
            #%% Inputting into arrays
            
            cmip6_deltaCs[model_i] = delta_Cs
            cmip6_deltaCstau[model_i] = deltaCs_tau
            cmip6_deltaCsRh[model_i] = deltaCs_rh
            cmip6_deltadelta[model_i] = deltadelta
            cmip6_remaining[model_i] = delta_Cs - deltaCs_tau - deltaCs_rh - deltadelta
            cmip6_NEPtau[model_i] = deltaCs_NEPtau
            cmip6_deltaNEPdeltatau[model_i] = deltaNEPdeltatau
            cmip6_NEPdeltatau[model_i] = deltaCs_NEPdeltatau

            cmip6_deltaCs_frac[model_i] = delta_Cs/soil_carbon_historical_data
            cmip6_deltaCstau_frac[model_i] = deltaCs_tau/soil_carbon_historical_data
            cmip6_deltaCsRh_frac[model_i] = deltaCs_rh/soil_carbon_historical_data
            cmip6_deltadelta_frac[model_i] = deltadelta/soil_carbon_historical_data
            
            
    #%% Saving data
    np.save("saved_data/cmip6_Cs_"+ssp+"_npp.npy", cmip6_deltaCs.data)
    np.save("saved_data/cmip6_Cstau_"+ssp+"_npp.npy", cmip6_deltaCstau.data)
    np.save("saved_data/cmip6_CsRh_"+ssp+"_npp.npy", cmip6_deltaCsRh.data)
    np.save("saved_data/cmip6_deltadelta_"+ssp+"_npp.npy", cmip6_deltadelta.data)
    np.save("saved_data/cmip6_remaining_"+ssp+"_npp.npy", cmip6_remaining.data)
    np.save("saved_data/cmip6_NEPtau_"+ssp+".npy", cmip6_NEPtau.data)
    np.save("saved_data/cmip6_deltaNEPdeltatau_"+ssp+".npy", cmip6_deltaNEPdeltatau.data)
    np.save("saved_data/cmip6_NEPdeltatau_"+ssp+".npy", cmip6_NEPdeltatau.data)
    
    np.save("saved_data/cmip6_fractionalCs_"+ssp+"_npp.npy", cmip6_deltaCs_frac.data)
    np.save("saved_data/cmip6_fractionalCstau_"+ssp+"_npp.npy", cmip6_deltaCstau_frac.data)
    np.save("saved_data/cmip6_fractionalCsRh_"+ssp+"_npp.npy", cmip6_deltaCsRh_frac.data)
    np.save("saved_data/cmip6_fractionaldeltadelta_"+ssp+"_npp.npy", cmip6_deltadelta_frac.data)