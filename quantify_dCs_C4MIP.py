#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script quantifies and saves dCs and the breakdown values (dCs,npp / dCs,tau / dnppdtau, dCs,nep / dCs,tau_nep / dnepdtau)
for the C4MIP simulations for 2xCO2 and 4xCO2.
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


#%%

# C4MIP simulations
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad']
C4MIP_labels = ['1pctCO2', '1pctCO2_bgc', '1pctCO2_rad']

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)

# 2xCO2 (65:70)
# 4xCO2 (135:140)


# Loop through each C4MIP simulation
for c4mip_option in range(0, len(C4MIP_simulation)):
    #
    c4mip = C4MIP_simulation[c4mip_option]
    C4MIP_label = C4MIP_labels[c4mip_option]
    
    
    # Empty arrays
    cmip6_deltaCs = np.zeros((len(cmip6_models)))
    cmip6_deltaCstau = np.zeros((len(cmip6_models)))
    cmip6_deltaCsNPP = np.zeros((len(cmip6_models)))
    cmip6_deltadelta = np.zeros((len(cmip6_models)))
    cmip6_NEPtau = np.zeros((len(cmip6_models)))
    cmip6_deltaNEPdeltatau = np.zeros((len(cmip6_models)))
    cmip6_NEPdeltatau = np.zeros((len(cmip6_models)))
    #
#    cmip6_deltaCs_fractional = np.zeros((len(cmip6_models)))
#    cmip6_deltaCstau_fractional = np.zeros((len(cmip6_models)))
#    cmip6_deltaCsNPP_fractional = np.zeros((len(cmip6_models)))
#    cmip6_deltadelta_fractional = np.zeros((len(cmip6_models)))
#    cmip6_NEPtau_fractional = np.zeros((len(cmip6_models)))
#    cmip6_deltaNEPdeltatau_fractional = np.zeros((len(cmip6_models)))
#    cmip6_NEPdeltatau_fractional = np.zeros((len(cmip6_models)))


    # for loop for each CMIP6 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i]
            print(model, c4mip)
            
            # land fraction
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)
            
            # Soil Carbon (cSoil)
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cSoil_Emon_'+model+'_'+c4mip+'*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = annual_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            
            # time dimension
            if model == 'ACCESS-ESM1-5':
                time_dimension = cSoil_cube.coord('year').points
        
            # Litter Carbon (cLitter)
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
            
            # NPP
            npp_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/npp_Lmon_'+model+'_'+c4mip+'*', model)
            npp_cube = open_netCDF(npp_cube)
            npp_cube = annual_average(npp_cube)
            npp_cube = global_total_percentage(npp_cube, landfrac=landfraction, latlon_cons=None)
            npp_data = npp_cube.data*86400.*360.            
            
            # tau 
            rh_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/rh_Lmon_'+model+'_'+c4mip+'*', model)
            rh_cube = open_netCDF(rh_cube)
            rh_cube = annual_average(rh_cube)
            rh_cube = global_total_percentage(rh_cube, landfrac=landfraction, latlon_cons=None)
            rh_data = rh_cube.data*86400.*360.
            #
            tau_data = Cs_data/rh_data

            # dCs
            deltaCs = np.mean(Cs_data[65:70]) - Cs_data[0]            
            
            # dCs,tau
            delta_tau = np.mean(tau_data[65:70]) - tau_data[0]
            deltaCs_tau = delta_tau*npp_data[0]
            
            # dCs,npp
            delta_npp = np.mean(npp_data[65:70]) - npp_data[0]
            deltaCs_npp = delta_npp*tau_data[0]
            
            # dnpp*dtau
            deltadelta = delta_npp*delta_tau
            
            # dNEP
            NEP_t = np.mean(npp_data[65:70]) - np.mean(rh_data[65:70])
            NEP_0 = npp_data[0] - rh_data[0]
            # dCs,NEP
            deltaCs_NEPtau = (NEP_t - NEP_0)*tau_data[0]
            # dCs,tau_NEP
            deltaCs_NEPdeltatau = NEP_0*delta_tau
            # dNEP*dtau
            deltaNEPdeltatau = (NEP_t - NEP_0)*delta_tau
            
            
            #%% Input to arrays
            cmip6_deltaCs[model_i] = deltaCs
            cmip6_deltaCstau[model_i] = deltaCs_tau
            cmip6_deltaCsNPP[model_i] = deltaCs_npp
            cmip6_deltadelta[model_i] = deltadelta
            cmip6_NEPtau[model_i] = deltaCs_NEPtau
            cmip6_NEPdeltatau[model_i] = deltaCs_NEPdeltatau
            cmip6_deltaNEPdeltatau[model_i] = deltaNEPdeltatau
            
#            cmip6_deltaCs_fractional[model_i] = deltaCs / Cs_data[0]
#            cmip6_deltaCstau_fractional[model_i] = deltaCs_tau / Cs_data[0]
#            cmip6_deltaCsNPP_fractional[model_i] = deltaCs_npp / Cs_data[0]
#            cmip6_deltadelta_fractional[model_i] = deltadelta / Cs_data[0]
#            cmip6_NEPtau_fractional[model_i] = deltaCs_NEPtau / Cs_data[0]
#            cmip6_NEPdeltatau_fractional[model_i] = deltaCs_NEPdeltatau / Cs_data[0]
#            cmip6_deltaNEPdeltatau_fractional[model_i] = deltaNEPdeltatau / Cs_data[0]


    #%% Saving data
    
    np.save('saved_data/cmip6_Cs_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCs.data)
    np.save('saved_data/cmip6_Cstau_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCstau.data)
    np.save('saved_data/cmip6_CsNPP_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCsNPP.data)
    np.save('saved_data/cmip6_deltadelta_'+C4MIP_label+'_2xCO2.npy', cmip6_deltadelta.data)
    np.save('saved_data/cmip6_NEPtau_'+C4MIP_label+'_2xCO2.npy', cmip6_NEPtau.data)
    np.save('saved_data/cmip6_deltaNEPdeltatau_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaNEPdeltatau.data)
    np.save('saved_data/cmip6_NEPdeltatau_'+C4MIP_label+'_2xCO2.npy', cmip6_NEPdeltatau.data)

#    np.save('saved_data/cmip6_fractionalCs_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCs_fractional.data)
#    np.save('saved_data/cmip6_fractionalCstau_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCstau_fractional.data)
#    np.save('saved_data/cmip6_fractionalCsNPP_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaCsNPP_fractional.data)
#    np.save('saved_data/cmip6_fractionaldeltadelta_'+C4MIP_label+'_2xCO2.npy', cmip6_deltadelta_fractional.data)
#    np.save('saved_data/cmip6_fractionalNEPtau_'+C4MIP_label+'_2xCO2.npy', cmip6_NEPtau_fractional.data)
#    np.save('saved_data/cmip6_fractionaldeltaNEPdeltatau_'+C4MIP_label+'_2xCO2.npy', cmip6_deltaNEPdeltatau_fractional.data)
#    np.save('saved_data/cmip6_fractionalNEPdeltatau_'+C4MIP_label+'_2xCO2.npy', cmip6_NEPdeltatau_fractional.data)
    
    