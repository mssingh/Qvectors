#
# Script to plot ERA5 data for EAE2122 assignment
#
# Choose the date and region below.
#

### Load the modules required ########################
print('Loading modules...')

# Plotting
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from matplotlib import colorbar
import matplotlib.patches as patches

# Meteorology functions from metpy
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT

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

# Function to make the map
from map_utils import make_map


### Input Parameters ########################################

# Time range to be plotted
t_start = "2025-10-25 06:00"
#t_end = "2025-10-25 06:00"
t_end = "2025-10-26 18:00"

timestamps = pd.date_range(start=t_start, end=t_end, freq="3h")


# Regions to be plotted in a lat-lon box
lats = [-45, -30]
lons = [135, 155]

level_q=850
level_sh=600

### Settings ################################################


# Variables we need
variables = ['cape','cin','u','v','q','t','u10','v10']


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
    #ds = ds.coarsen(latitude=4, longitude=4, boundary='trim').mean()


    ### Calculate variables #######################################
    print('    Calculating...')


    plt.close('all')


    du = ds.u.sel(level=level_sh)-ds.u.sel(level=1000)
    dv = ds.v.sel(level=level_sh)-ds.v.sel(level=1000)

    shear = (du**2 + dv**2)**0.5

    u = ds.u10
    v = ds.v10
    div = mpcalc.divergence(u,v)
 
    # Smooth the divergence
    div = mpcalc.smooth_gaussian(div, 8)
    #div = mpcalc.smooth_n_point(div, 9,5) 

    q = ds.q.sel(level=level_q)
    t = ds.t.sel(level=level_q)

    rh = mpcalc.relative_humidity_from_specific_humidity(level_q*units.hPa,t,q)

    cape = ds.cape
    cin = ds.cin


    ### Make the plots ###########################################
    print('    PLotting...')

    ## Plot the CAPE and CIN

    # Make the map
    fig,ax = make_map(lons,lats,grid_spacing=5,states=True,Melbourne=True,regional=True)


    # Plot CAPE
    cmap = plt.cm.Reds
    clevs = np.arange(0, 1601, 200)
    cf = ax.contourf(ds.longitude, ds.latitude, cape, clevs, cmap=cmap,
                 extend='max')

    # Make colorbar
    cb = plt.colorbar(cf, location='right', orientation="vertical",pad=0.05, aspect=50,ticks=clevs,shrink=0.5)
    cb.ax.tick_params(labelsize=14)
    cb.set_ticks([0, 200, 400, 600, 800, 1000,1200,1400,1600])
    cb.set_label('CAPE (J/kg)',fontsize=14)


    # Plot CIN
    cmap = plt.cm.winter
    clevs = np.arange(100, 501, 100)
    cs = ax.contour(ds.longitude, ds.latitude, cin, clevs, cmap=cmap,extend='both')

    # Make colorbar
    norm = colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # required even though unused

    cb2 = plt.colorbar(sm, ax=ax,location='left', orientation="vertical", pad=0.08, aspect=50,ticks=clevs,shrink=0.5,extend="both")
    cb2.ax.tick_params(labelsize=14)
    cb2.set_ticks([100, 200,300,400,500])
    cb2.set_label('CIN (J/kg)',fontsize=14)


    # Add some titles
    plt.title('CAPE (shading; J/kg), CIN (contours; J/kg)', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_CAPE_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(fig)



    ## Plot the humidity and shear

    fig,ax = make_map(lons,lats,grid_spacing=5,states=True,Melbourne=True,regional=True)


    # Plot humdity
    cmap = plt.cm.Blues
    clevs = np.arange(20, 130, 10)
    cf = ax.contourf(ds.longitude, ds.latitude, rh*100, clevs, cmap=cmap,extend='both')

    # Make colorbar
    cb = plt.colorbar(cf, location='right', orientation="vertical",pad=0.05, aspect=50,ticks=clevs,shrink=0.5)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('relative humidity (%)',fontsize=14)
    cb.ax.set_ylim(20, 100)


    # Plot shear
    cmap = plt.cm.plasma
    clevs = np.arange(0, 40, 5)
    cs = ax.contour(ds.longitude, ds.latitude, shear, clevs, cmap=cmap)
    ax.clabel(cs, fmt='%d',fontsize=14)

    # Make colorbar
    norm = colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # required even though unused

    cb2 = plt.colorbar(sm, ax=ax,location='left', orientation="vertical", pad=0.08, aspect=50,ticks=clevs,shrink=0.5,extend="both")
    cb2.ax.tick_params(labelsize=14)
    cb2.set_label('0-6 km shear (m/s)',fontsize=14)


    # Add some titles
    plt.title('850 hPa relative humidity (shading; %), 0-6km shear (contours)', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_humidity_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(fig)



    ## Plot the winds and divergence

    fig,ax = make_map(lons,lats,grid_spacing=5,states=True,Melbourne=True,regional=True)


    # Plot divergence
    cmap = plt.cm.coolwarm
    clevs = np.arange(-10,10,0.5)
    cf = ax.contourf(ds.longitude, ds.latitude, div*100000, clevs, cmap=cmap,extend='both')

    # Make colorbar
    cb = plt.colorbar(cf, location='right', orientation="vertical",pad=0.05, aspect=50,ticks=clevs,shrink=0.5)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('divergnce ($10^{-5}$ s$^{-1}$)',fontsize=14)
    cb.set_ticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])


    # Plot surf winds
    cf2 = ax.quiver(ds.longitude[::2], ds.latitude[::2], u[::2,::2], v[::2,::2],scale=250)

    cb2 = plt.colorbar(cf, ax=ax,location='left', orientation="vertical", pad=0.08, aspect=50,ticks=clevs,shrink=0.5,extend="both")
    cb2.ax.tick_params(labelsize=14)
    cb2.set_label('divergnce ($10^{-5}$ s$^{-1}$)',fontsize=14)
    cb2.set_ticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])
    #cb.axis('off')


    # Add some titles
    plt.title('near-surface winds (arrows) and divergence (colours; $10^{-5}$ s$^{-1}$)', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    cb2.ax.set_visible(False)

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_surf-winds_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    #plt.show()
    plt.close(fig)



    ## Plot the soundings

    # Get the melbroune profile
    p = ds.level
    T = ds.t.sel(latitude=-37.81,longitude=144.96,method='nearest')-273.15
    Q = ds.q.sel(latitude=-37.81,longitude=144.96,method='nearest')
    Td = mpcalc.dewpoint_from_specific_humidity(p, Q) 
    u = ds.u.sel(latitude=-37.81,longitude=144.96,method='nearest')
    v = ds.v.sel(latitude=-37.81,longitude=144.96,method='nearest')


    # Set up the Skew-T
    fig = plt.figure(figsize=(12, 12))
    skew = SkewT(fig, rotation=45)

    skew.plot(p, T, 'r', linewidth=2)
    skew.plot(p, Td, 'g', linewidth=2)
    skew.plot_barbs(p, u, v)

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 60)
    skew.ax.set_xlabel('Temperature (°C)',fontsize=14)
    skew.ax.set_ylabel('Pressure (hPa)',fontsize=14)
    skew.ax.tick_params(axis='both', labelsize=14)  

    # Fiducial lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    # Reverse the martices
    p = p[::-1] 
    T = T[::-1] 
    Td = Td[::-1] 

    parcel_prof = mpcalc.parcel_profile(p, T[0].values*units("degC"), Td[0].values*units("degC"))

    parcel_prof = (parcel_prof.values-273.15)*units("degC")
    skew.plot(p, parcel_prof, 'k', linewidth=2)
    #skew.shade_cape(p, T, parcel_prof)
    #skew.shade_cin(p, T, parcel_prof, Td)
 


    # Add some titles
    plt.title('Melbourne forecast sounding', loc='left',fontsize=16)
   
    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()

    time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/ERA5_skewT_"+time_str+".png"

    plt.savefig(fname, dpi=200, bbox_inches="tight")
    #plt.show()
    plt.close(fig)





