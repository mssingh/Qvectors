#
# Script to plot ERA5 data for EAE2122 assignment
#
# Choose the date and region below.
#

### Load the modules required ########################
print('Loading modules...')

# PLotting and Mapping
import cartopy.crs as ccrs          
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import colors


# Meteorology functions from metpy
import metpy.calc as mpcalc
from metpy.units import units

# Numerics and data
import numpy as np
import xarray as xr
import pandas as pd

# General utilities
import os
import glob
from datetime import datetime
import subprocess
from pathlib import Path
import shutil

import warnings
warnings.filterwarnings("ignore", message="facecolor will have no effect")

# Function to read the ERA5 data into a dataset 
from era5_utils import open_era5_month
from map_utils import make_map


### Input Parameters ########################################

# Time range to be plotted
t_start = "2025-10-25 00:00"
#t_end = "2025-10-25 00:00"
t_end = "2025-10-27 00:00"

timestamps = pd.date_range(start=t_start, end=t_end, freq="6h")


# Regions to be plotted in a lat-lon box
lats = [-55, -10]
lons = [90, 175]

level_t=850
level_u=200

### Settings ################################################


# Variables we need
variables = ['z','u','v','t','msl','tp']



# Set subset slice for the geographic extent of data to limit download
lon_slice = slice(lons[0]-4,lons[1]+4)
lat_slice = slice(lats[0]-4,lats[1]+4)

### Loop over each time ####################################
print('Looping over times...')


i = 0
for tt in timestamps:
    i = i+1
    print(tt)

    ### Load data ###############################################
    print('    Loading data...')

    # Get the variables we require for the timeslice
    ds = open_era5_month(tt,variables=variables)
    ds = ds.squeeze()



    # Ensure the data is sorted with increasing lon and lat
    ds = ds.sortby('latitude')
    ds = ds.sortby('longitude')
   
    ds = ds.sel(latitude=lat_slice, longitude=lon_slice)


    # Coarsen to 1x1 deg
    ds = ds.coarsen(latitude=4, longitude=4, boundary='trim').mean()

    ### Calculate variables #######################################
    print('    Calculating...')


    u = ds.u.sel(level=level_u)
    v = ds.v.sel(level=level_u)
    z = ds.z.sel(level=level_u)/9.81
    t = ds.t.sel(level=level_t)

    dz = (ds.z.sel(level=500)-ds.z.sel(level=1000))/9.81

    speed = (u**2 + v**2)**0.5

    msl = ds.msl

    ### Make the plots ##############################################
    print('    Plotting...')

    ## Plot temperature and MSLP

    # Make the map
    fig,ax = make_map(lons,lats,grid_spacing=10,states=False,Melbourne=False,W=14,H=12,regional=True)


    # Plot 850-hPa Temperatures
    cmap = plt.cm.turbo
    clevs_t = np.arange(250, 300, 2)
    cf = ax.contourf(ds.longitude, ds.latitude, t, clevs_t, cmap=cmap, extend='both')

    # Make a colorbar
    cb = plt.colorbar(cf, orientation='horizontal', pad=0.075, aspect=50,ticks=clevs_t,shrink=0.7)
    cb.ax.tick_params(labelsize=14)
    cb.set_ticks([250, 260, 270, 280, 290, 300])
    cb.set_label('850 hPa temperature (K)',fontsize=14)


    # Plot Mean sea-level pressure
    clevs_msl = np.arange(920, 1040, 4)
    cs = ax.contour(ds.longitude, ds.latitude, msl/100, clevs_msl, colors='black')
    plt.clabel(cs, fmt='%d')


    # Add some titles
    plt.title('MSLP (black; hPa), 850 hPa temperature (colours; K)', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_MSLP_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    #plt.show()
    plt.close(fig)

    ## Plot windspeed and geopotenial at 200 hPa

    # Make the map
    fig,ax = make_map(lons,lats,grid_spacing=10,states=False,Melbourne=False,W=14,H=12,regional=True)



    # Plot 200-hPa wind speed
    cmap = plt.cm.Greens
    clevs_t = np.arange(40, 75, 5)
    cf = ax.contourf(ds.longitude, ds.latitude, speed, clevs_t, cmap=cmap,extend='both')

    # Make colorbar
    cb = plt.colorbar(cf, orientation='horizontal', pad=0.075, aspect=50,ticks=clevs_t,shrink=0.7)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('200 hPa windspeed (m/s)',fontsize=14)


    # Plot 200 hPa geopotential
    clevs_z = np.arange(10, 15, 0.1)
    cs = ax.contour(ds.longitude, ds.latitude, z/1000, clevs_z, colors='black')

 
    plt.clabel(cs)#levels=[80,90,100,110,120,130,140,150])


    # Add some titles
    plt.title('200 hPa Geopotential height (black; km), 200 hPa windspeed (colours)', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_200hPa_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    #plt.show()
    plt.close(fig)







