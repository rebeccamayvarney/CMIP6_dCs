#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 2022
@author: Rebecca Varney, University of Exeter (r.varney@exeter.ac.uk)

"
Script plots dCs maps for each CMIP6 ESMs in SSP585 simulations.
"

"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, open_netCDF, select_time, time_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%% Set up subplot figure

fig = plt.figure(1, figsize=(46,52))
gs = gspec.GridSpec(5, 3, figure=fig, width_ratios=[1, 1, 0.075], hspace=0.4, wspace=0.3)
column_1 = 0
row_1 = 0
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
    'axes.labelsize':58,
    'xtick.labelsize':58,
    'ytick.labelsize':58,
    'font.size':58,
}
plt.rcParams.update(params)


#%%

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)

# SSP scenario
ssp='ssp585' # 'ssp126', 'ssp245'


# for loop for each CMIP6 model
for model_i in range(0, n_models):
    model = cmip6_models[model_i] # seleting the models
    print(ssp, model)
 
    if ssp=='ssp245' or ssp=='ssp585':
        if model=='GFDL-ESM4':
            continue

    
    # cSoil (historical)
    cSoil_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_historical*', model)
    cSoil_historical_cube = open_netCDF(cSoil_historical_cube)
    n_lat = cSoil_historical_cube.coord('latitude').points
    n_lon = cSoil_historical_cube.coord('longitude').points
    cSoil_historical_cube = select_time(cSoil_historical_cube, 1995, 2005)   
    cSoil_historical_cube = time_average(cSoil_historical_cube)
    cSoil_historical_data = cSoil_historical_cube.data
    
    # cSoil (future)
    cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cSoil_Emon_'+model+'_'+ssp+'_*', model)
    cSoil_cube = open_netCDF(cSoil_cube)
    cSoil_cube = select_time(cSoil_cube, 2090, 2100)
    cSoil_cube = time_average(cSoil_cube)
    cSoil_data = cSoil_cube.data


    if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
        # cLitter (historical)
        cLitter_historical_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_historical*', model)
        cLitter_historical_cube = open_netCDF(cLitter_historical_cube)     
        cLitter_historical_cube = select_time(cLitter_historical_cube, 1995, 2005)   
        cLitter_historical_cube = time_average(cLitter_historical_cube)
        cLitter_historical_data = cLitter_historical_cube.data
        soil_carbon_historical_data = cSoil_historical_data + cLitter_historical_data
        # cLitter (future)
        cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/cmip6_data/cLitter_Lmon_'+model+'_'+ssp+'*', model)
        cLitter_cube = open_netCDF(cLitter_cube)
        cLitter_cube = select_time(cLitter_cube, 2090, 2100)
        cLitter_cube = time_average(cLitter_cube)
        cLitter_data = cLitter_cube.data
        #
        soil_carbon_future_data = cSoil_data + cLitter_data
    else:
        soil_carbon_historical_data = cSoil_historical_data.copy()
        soil_carbon_future_data = cSoil_data.copy()
    
    # dCs
    delta_Cs = soil_carbon_future_data - soil_carbon_historical_data
    delta_Cs = ma.masked_where(delta_Cs==0, delta_Cs)


    #%% Plotting
    print(row_1, column_1)
    ax = fig.add_subplot(gs[row_1, column_1], projection=ccrs.PlateCarree())
    
    #  add lat/lon grid lines to the figure
    gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
    gl.yformatter=LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60])
    gl.xlocator = mticker.FixedLocator([-90, 0, 90])
    gl.ylabels_left=True
    gl.xformatter=LONGITUDE_FORMATTER
    gl.xlabels_bottom=True
    ax.coastlines()
    # set up the x and y coordination
    lat = n_lat
    lon = n_lon
    x, y = np.meshgrid(lon, lat)
    
    #print(np.min(delta_Cs), np.max(delta_Cs))
    line = np.arange(-10, 10, 0.5)
    diff = plt.contourf(x, y, delta_Cs, line, cmap='BrBG_r', extend='both', transform=ccrs.PlateCarree(central_longitude=0))
    ax.set_title(model)
    ax.set_ylim(-70,90)
    
    
    #%% increase row and column 
    row_1 += 1 
    if row_1==5:
        column_1 += 1
        row_1 = 0


#%%

# Colour bar
ax=fig.add_subplot(gs[:,2])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$\Delta C_{s}$ (kg C m$^{-2}$)')

# Save figure
fig.savefig('figures/dCs_maps_cmip6_SSP585_v1', bbox_inches='tight')
plt.close()