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


### Input Parameters ########################################

# Time range to be plotted
t_start = "2025-10-25 00:00"
t_end = "2025-10-27 00:00"

timestamps = pd.date_range(start=t_start, end=t_end, freq="6h")


# Regions to be plotted in a lat-lon box
lats = [-10, -55]
lons = [90, 175]

level_t=850
level_u=200

### Settings ################################################


# Variables we need
variables = ['z','u','v','t','msl','tp']



# Set subset slice for the geographic extent of data to limit download
lon_slice = slice(lons[0]-4,lons[1]+4)
lat_slice = slice(lats[0]+4,lats[1]-4)

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
   

    # Coarsen to 1x1 deg
    ds = ds.coarsen(latitude=4, longitude=4, boundary='trim').mean()


    ### Make the plots ##############################################
    print('    Plotting...')


    u = ds.u.sel(level=level_u)
    v = ds.v.sel(level=level_u)
    z = ds.z.sel(level=level_u)/9.81
    t = ds.t.sel(level=level_t)

    dz = (ds.z.sel(level=500)-ds.z.sel(level=1000))/9.81

    speed = (u**2 + v**2)**0.5

    msl = ds.msl


    ## Plot number 1


    # Set the map projection (how the data will be displayed)
    mapcrs = ccrs.PlateCarree()

    # Set the data projection (GFS is lat/lon format)
    datacrs = ccrs.PlateCarree()

    # Start the figure and set an extent to only display a smaller graphics area
    fig = plt.figure(1, figsize=(14, 12))
    ax = plt.subplot(111, projection=mapcrs)
    ax.set_extent([lons[0], lons[1], lats[1], lats[0]], ccrs.PlateCarree())

    # Add map features to plot coastlines and state boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),color='black')

    # Add gridlines with nicely spaced labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 14}
    gl.ylabel_style = {"size": 14}
    gl.xlocator = plt.FixedLocator(range(int(np.floor(lons[0]/10)*10),int( np.ceil(lons[1]/10)*10), 10))  # every 10 degrees longitude
    gl.ylocator = plt.FixedLocator(range(int(np.ceil(lats[1]/10)*10), int(np.floor(lats[0]/10)*10), 10))   # every 10 degrees latitude


    # Create colormap: 
    #cmap = plt.cm.gist_ncar
    cmap = plt.cm.turbo
    #cmap = plt.cm.viridis

   # Plot 850-hPa Temperatures
    clevs_t = np.arange(250, 300, 2)
    cf = ax.contourf(ds.longitude, ds.latitude, t, clevs_t, cmap=cmap,
                 extend='both', transform=datacrs)
    cb = plt.colorbar(cf, orientation='horizontal', pad=0.075, aspect=50,
                  ticks=clevs_t,shrink=0.7)
    cb.ax.tick_params(labelsize=14)
    cb.set_ticks([250, 260, 270, 280, 290, 300])
    cb.set_label('850 hPa temperature (K)',fontsize=14)


    # Plot Mean sea-level pressure
    clevs_msl = np.arange(920, 1040, 4)
    cs = ax.contour(ds.longitude, ds.latitude, msl/100, clevs_msl, colors='black', transform=datacrs)
    plt.clabel(cs, fmt='%d')

    # PLot thickness
    #clevs_dz = np.arange(400,600,4)
    #csf = ax.contour(ds.longitude, ds.latitude, dz/10, clevs_dz, colors='black', transform=datacrs)
    #plt.clabel(csf, fmt='%d')



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

    ## Plot number 2


    # Set the map projection (how the data will be displayed)
    mapcrs = ccrs.PlateCarree()

    # Set the data projection (GFS is lat/lon format)
    datacrs = ccrs.PlateCarree()

    # Start the figure and set an extent to only display a smaller graphics area
    fig = plt.figure(1, figsize=(14, 12))
    ax = plt.subplot(111, projection=mapcrs)
    ax.set_extent([lons[0], lons[1], lats[1], lats[0]], ccrs.PlateCarree())

    # Add map features to plot coastlines and state boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),color='black')

    # Add gridlines with nicely spaced labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 12}
    gl.ylabel_style = {"size": 12}
    gl.xlocator = plt.FixedLocator(range(int(np.floor(lons[0]/10)*10),int( np.ceil(lons[1]/10)*10), 10))  # every 10 degrees longitude
    gl.ylocator = plt.FixedLocator(range(int(np.ceil(lats[1]/10)*10), int(np.floor(lats[0]/10)*10), 10))   # every 10 degrees latitude


    # Create colormap: 
    cmap = plt.cm.Greens

   # Plot 200-hPa wind speed
    clevs_t = np.arange(40, 75, 5)
    cf = ax.contourf(ds.longitude, ds.latitude, speed, clevs_t, cmap=cmap,
                 extend='both', transform=datacrs)
    cb = plt.colorbar(cf, orientation='horizontal', pad=0.075, aspect=50,
                  ticks=clevs_t,shrink=0.7)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('200 hPa windspeed (m/s)',fontsize=14)

   # Plot 200-hPa winds
   # clevs_t = np.arange(250, 300, 2)
   # cf = ax.quiver(ds.longitude, ds.latitude, u, v, transform=datacrs)


    # Plot 200 hPa geopotential
    clevs_z = np.arange(10, 15, 0.1)
    cs = ax.contour(ds.longitude, ds.latitude, z/1000, clevs_z, colors='black', transform=datacrs)

 
    plt.clabel(cs)#levels=[80,90,100,110,120,130,140,150])

    # PLot thickness
    #clevs_dz = np.arange(400,600,4)
    #csf = ax.contour(ds.longitude, ds.latitude, dz/10, clevs_dz, colors='black', transform=datacrs)
    #plt.clabel(csf, fmt='%d')



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







