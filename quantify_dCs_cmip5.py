#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 2021
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script quantifies and saves dCs and the breakdown values (dCs,npp / dCs,tau / dnppdtau, dCs,nep / dCs,tau_nep / dnepdtau)
for CMIP5 ESMs for RCP simulations.
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
from rmv_cmip_analysis import combine_netCDF_cmip5, combine_netCDF_model, open_netCDF
from rmv_cmip_analysis import select_time, time_average, global_total_percentage


#%% CMIP5

# rcp scenarios
rcp_options = ['rcp26', 'rcp45', 'rcp85']
rcp_options_length = len(rcp_options)

# CMIP5 models
cmip5_models = ['BNU-ESM', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)


#%%
# Loop through each rcp run being considered
for rcp_option in range(0, rcp_options_length):
    rcp = rcp_options[rcp_option] # selecting the rcp scenario

    # Empty arrays
    cmip5_deltaCs = np.zeros((len(cmip5_models)))
    cmip5_deltaCstau = np.zeros((len(cmip5_models)))
    cmip5_deltaCsRh = np.zeros((len(cmip5_models)))
    cmip5_deltadelta = np.zeros((len(cmip5_models)))
    cmip5_remaining = np.zeros((len(cmip5_models)))
    cmip5_percent = np.zeros((len(cmip5_models)))
    cmip5_NEPtau = np.zeros((len(cmip5_models)))
    cmip5_deltaNEPdeltatau = np.zeros((len(cmip5_models)))
    cmip5_NEPdeltatau = np.zeros((len(cmip5_models)))
    #
    cmip5_deltaCs_frac = np.zeros((len(cmip5_models)))
    cmip5_deltaCstau_frac = np.zeros((len(cmip5_models)))
    cmip5_deltaCsRh_frac = np.zeros((len(cmip5_models)))
    cmip5_deltadelta_frac = np.zeros((len(cmip5_models)))


    # for loop for each CMIP5 model
    for model_i in range(0, n_models):
            model = cmip5_models[model_i] # seleting the models
            print(rcp, model)
            
            # land fraction
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip5_data/sftlf_fx_'+model+'_*', model)
            
            
            #%%
            # cSoil (historical)
            cSoil_historical_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cSoil_Lmon_'+model+'_historical*', 'soil_carbon_content', model)
            cSoil_historical_cube = open_netCDF(cSoil_historical_cube)
            cSoil_historical_cube = select_time(cSoil_historical_cube, 1995, 2005)
            cSoil_historical_cube = time_average(cSoil_historical_cube)
            cSoil_historical_cube = global_total_percentage(cSoil_historical_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_historical_data = cSoil_historical_cube.data
            
            # cLitter (historical)
            if model=='BNU-ESM' or model=='CESM1-CAM5' or model=='CanESM2' or model=='IPSL-CM5A-LR' or model=='MIROC-ESM' or model=='MPI-ESM-LR' or model=='NorESM1-M':
                cLitter_historical_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/cLitter_Lmon_'+model+'_historical*', 'litter_carbon_content', model)
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
            rh_historical_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/rh_Lmon_'+model+'_historical*', 'heterotrophic_respiration_carbon_flux', model)
            rh_historical_cube = open_netCDF(rh_historical_cube)
            rh_historical_cube = select_time(rh_historical_cube, 1995, 2005)
            rh_historical_cube = time_average(rh_historical_cube)
            rh_historical_cube = global_total_percentage(rh_historical_cube, landfrac=landfraction, latlon_cons=None)
            rh_historical_data = rh_historical_cube.data*86400.*360.
            #
            tau_0_data = soil_carbon_historical_data/rh_historical_data
            
            # NPP (historical)
            npp_historical_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_data/npp_Lmon_'+model+'_historical*', 'net_primary_productivity_of_carbon', model)
            npp_historical_cube = open_netCDF(npp_historical_cube)
            npp_historical_cube = select_time(npp_historical_cube, 1995, 2005)
            npp_historical_cube = time_average(npp_historical_cube)
            npp_historical_cube = global_total_percentage(npp_historical_cube, landfrac=landfraction, latlon_cons=None)
            npp_historical_data = npp_historical_cube.data*86400.*360.
            
            
            # cSoil (future)
            cSoil_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_historical_'+rcp+'/cSoil_Lmon_'+model+'_*', 'soil_carbon_content', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = select_time(cSoil_cube, 2090, 2100)
            cSoil_cube = time_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            
            # cLitter (future)
            if model=='BNU-ESM' or model=='CESM1-CAM5' or model=='CanESM2' or model=='IPSL-CM5A-LR' or model=='MIROC-ESM' or model=='MPI-ESM-LR' or model=='NorESM1-M':
                cLitter_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_historical_'+rcp+'/cLitter_Lmon_'+model+'_*', 'litter_carbon_content', model)
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
            rh_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_historical_'+rcp+'/rh_Lmon_'+model+'_*', 'heterotrophic_respiration_carbon_flux', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = select_time(rh_cube, 2090, 2100)
            rh_cube = time_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube, landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data*86400.*360.
            #
            tau_data = soil_carbon_future_data/rh_data
            
            # NPP (future)
            npp_cube = combine_netCDF_cmip5('/home/rmv203/DATA/cmip5_historical_'+rcp+'/npp_Lmon_'+model+'_*', 'net_primary_productivity_of_carbon', model)
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
            
            cmip5_deltaCs[model_i] = delta_Cs
            cmip5_deltaCstau[model_i] = deltaCs_tau
            cmip5_deltaCsRh[model_i] = deltaCs_rh
            cmip5_deltadelta[model_i] = deltadelta
            cmip5_remaining[model_i] = delta_Cs - deltaCs_tau - deltaCs_rh - deltadelta
            cmip5_NEPtau[model_i] = deltaCs_NEPtau
            cmip5_deltaNEPdeltatau[model_i] = deltaNEPdeltatau
            cmip5_NEPdeltatau[model_i] = deltaCs_NEPdeltatau
            
            cmip5_deltaCs_frac[model_i] = delta_Cs/soil_carbon_historical_data
            cmip5_deltaCstau_frac[model_i] = deltaCs_tau/soil_carbon_historical_data
            cmip5_deltaCsRh_frac[model_i] = deltaCs_rh/soil_carbon_historical_data
            cmip5_deltadelta_frac[model_i] = deltadelta/soil_carbon_historical_data
            
            
    #%% Saving data

    np.save("saved_data/cmip5_Cs_"+rcp+"_npp.npy", cmip5_deltaCs.data)
    np.save("saved_data/cmip5_Cstau_"+rcp+"_npp.npy", cmip5_deltaCstau.data)
    np.save("saved_data/cmip5_CsRh_"+rcp+"_npp.npy", cmip5_deltaCsRh.data)
    np.save("saved_data/cmip5_deltadelta_"+rcp+"_npp.npy", cmip5_deltadelta.data)
    np.save("saved_data/cmip5_remaining_"+rcp+"_npp.npy", cmip5_remaining.data)
    np.save("saved_data/cmip5_NEPtau_"+rcp+"_npp.npy", cmip5_NEPtau.data)
    np.save("saved_data/cmip5_deltaNEPdeltatau_"+rcp+".npy", cmip5_deltaNEPdeltatau.data)
    np.save("saved_data/cmip5_NEPdeltatau_"+rcp+".npy", cmip5_NEPdeltatau.data)
    
    np.save("saved_data/cmip5_fractionalCs_"+rcp+"_npp.npy", cmip5_deltaCs_frac.data)
    np.save("saved_data/cmip5_fractionalCstau_"+rcp+"_npp.npy", cmip5_deltaCstau_frac.data)
    np.save("saved_data/cmip5_fractionalCsRh_"+rcp+"_npp.npy", cmip5_deltaCsRh_frac.data)
    np.save("saved_data/cmip5_fractionaldeltadelta_"+rcp+"_npp.npy", cmip5_deltadelta_frac.data)