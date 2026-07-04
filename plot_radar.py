#
# Script to plot radar data for EAE2122 assignment
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

# Distance function
from pyproj import Geod

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


# Function to make the map
from map_utils import make_map


### Input Parameters ########################################

# Time range to be plotted
t_start = "2025-10-25 18:00"
t_end = "2025-10-26 18:00"
#t_end = "2025-10-26 18:00"

timestamps = pd.date_range(start=t_start, end=t_end, freq="5min")


# Radar ID
radar_ID = "2" # Melbourne
radar_lon = 144.7554
radar_lat = -37.8553



### Settings ################################################


### Loop over each time ####################################
print('Looping over times...')


i = 0
for tt in timestamps:
    i = i+1

    dt = pd.to_datetime(tt, unit='s', utc=True)
    time_str = dt.strftime("%Y-%m-%d_%H:%M")
    fname = "Figures/radar/radar_"+time_str+".png"

    if os.path.exists(fname):
        continue  # skip this iteration

    print(dt.strftime('%Y-%m-%d %H:%M:%S'))

    ### Load data ###############################################
    print('    Loading data...')

    # Get the radar snapshot
    date_string = dt.strftime('%Y%m%d_%H%M%S')
    filename = './radar/'+radar_ID+'_'+date_string+'.prcp-crate.nc'

    if not os.path.exists(filename):
        print("    Cannot find file... skipping")
        continue  # skip this iteration

    ds = xr.open_dataset(filename)    



    ### Calculate variables #######################################
    print('    Calculating...')

    # Calculate the latitude and longitude
    # This is relatively accurate, but ignores the difference in the longitude distance between the northern and southern edge of the radar
    geod = Geod(ellps='WGS84')
    n = len(ds.x)
    lon,yy,zz = geod.fwd(np.full(n,radar_lon), np.full(n,radar_lat), np.full(n,90), ds.x*1000)
    xx,lat,zz = geod.fwd(np.full(n,radar_lon), np.full(n,radar_lat), np.full(n,0), ds.y*1000)

    pr = ds.rain_rate

    # Limit to 128 km from radar
    X, Y = np.meshgrid(ds.x, ds.y)
    r = (X**2 + Y**2)**0.5
    #pr = pr.where(r<127.) 

    # Regions to be plotted in a lat-lon box
    lats = [min(lat), max(lat)]
    lons = [min(lon),max(lon)]

    ### Make the plots ###########################################
    print('    PLotting...')


    # Make the map
    fig,ax,mapcrs = make_map(lons,lats,grid_spacing=20,states=False,Melbourne='black',regional=True,hires=True)


    # Add some guiding lines
    clevs = [50,100]
    cs = ax.contour(lon,lat, r, clevs, colors="brown")
    ax.clabel(cs, fmt='%d',fontsize=14)

    ax.plot([radar_lon, radar_lon],[min(lat), max(lat)],color="brown",linewidth=1,transform=mapcrs)
    ax.plot([min(lon),max(lon)],[radar_lat, radar_lat],color="brown",linewidth=1,transform=mapcrs)


    # Add some towns
    suburbs = [
     "CBD",
     "Geelong",
     "Werribee",
     "Sunbury",
     "Dandenong",
     "Frankston",
     "Mornington",
     "Wallan",
     "Box Hill",
     "Preston",
     "Pakenham",
     "Melton",
     "Melbourne Ap",
     "Monash Uni",
     "Macedon", 
     "Healesville"
    ]

    lonp = [
    144.9631, 144.3547, 144.6570, 144.7280,
    145.1900, 145.1220, 145.0380,144.9840,
    145.1210, 145.0000, 145.4860, 144.5830,
    144.8410,145.1340,144.5620, 145.5260
    ]

    latp = [
    -37.8136, -38.1499, -37.9000, -37.5790, 
    -37.9830, -38.1446, -38.2175,-37.4160,
    -37.8190, -37.7400, -38.0770, -37.6830,
    -37.6690,-37.9078,-37.4230, -37.6530
    ]



    import cartopy.crs as ccrs
    for c, lo, la in zip(suburbs, lonp, latp):
       ax.scatter(lo,la, marker='o', color='gray', s=30,transform=ccrs.PlateCarree(),zorder=5)
       ax.text(lo+0.01,la+0.01,c,fontsize=10,color="black",transform=ccrs.PlateCarree())


    colorlist = [
    "#FFFFFF", "#F0FBFF", "#D6F3FF", "#BDEBFF", "#9FDFFF",
    "#7CCFFF", "#5BB8FF", "#3A9BFF", "#1F7DFF", "#0B5CFF",
    "#00E676", "#00D45A", "#00C853", "#00B84D", "#00A843",
    "#FFFF66", "#FFEB3B", "#FFD54F", "#FFC107", "#FFB300",
    "#FF9800", "#FF7A00", "#FF5A1F", "#FF3D2E", "#E53935",
    "#C62828", "#AD1457", "#8E24AA", "#5E35B1", "#000000"
    ]

    cmap = colors.ListedColormap(colorlist)



    # Plot rain_rate
    pr = pr.where(pr>0.5)
    clevs = np.arange(0, 30, 1)
    cf = ax.contourf(lon,lat, pr, clevs, cmap=cmap,
                 extend='max')

    # Make colorbar
    cb = plt.colorbar(cf, location='bottom', orientation="horizontal",pad=0.075, aspect=50,ticks=clevs,shrink=0.5)
    cb.ax.tick_params(labelsize=14)
    cb.set_ticks([0,5,10,15,20,25,30])
    cb.set_label('rainfall rate (mm/hr)',fontsize=14)

   


    # Add some titles
    plt.title('Melbourne radar', loc='left',fontsize=16)
   
    time_str = dt.strftime("%Y-%m-%d %H:%M")
    plt.title('{} UTC'.format(time_str), loc='right',fontsize=16)

    plt.tight_layout()


    plt.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(fig)






