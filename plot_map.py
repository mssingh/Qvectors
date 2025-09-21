#
# Script to plot QG vectors from ERA5 data
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
import glob
from datetime import datetime

# Function to read the ERA5 data into a dataset 
from era5_utils import open_era5_month



### Input Parameters ########################################

# Time of snapshot to be plotted
yyyy = 2022
mm = 9
dd = 24
hh = 0

# Regions to be plotted in a lat-lon box
lats = [-10, -55]
lons = [90, 175]


# Pressure level to be plotted
level = 850


### Settings ################################################


# Variables we need
variables = ['z','u','v','w','t','msl','pt']


# Smoothing for the calculation of Q-vectors 
# Note: we first coarsen the data to 1x1 deg
n_smooth = 9                # number of points of the smoother
N_smooth_q = 10             # numer of times applied for qvectors


### Load data ###############################################
print('Loading data...')

# Get the variables we require for the timeslice
ds = open_era5_month(yyyy, mm, day=dd,hour=hh,variables=variables)


# Set subset slice for the geographic extent of data to limit download
lon_slice = slice(lons[0]-4,lons[1]+4)
lat_slice = slice(lats[0]+4,lats[1]-4)

ds = ds.sel(longitude=lon_slice,latitude=lat_slice,level=level)

# Ensure the data is sorted with increasing lon and lat
ds = ds.sortby('latitude')
ds = ds.sortby('longitude')

# Coarsen to 1x1 deg
ds_orig = ds;
ds = ds.coarsen(latitude=4, longitude=4, boundary='trim').mean()


# Smooth variables used for Q-vector calculation
u = mpcalc.smooth_n_point(ds.u, n_smooth, N_smooth_q)
v = mpcalc.smooth_n_point(ds.v, n_smooth, N_smooth_q)
t = mpcalc.smooth_n_point(ds.t, n_smooth, N_smooth_q)



### Compute the Q-vectors #######################################
print('Computing Q-vectors...')

# Compute grid spacings for data
dx, dy = mpcalc.lat_lon_grid_deltas(ds.longitude,ds.latitude)

# Compute the Q-vector components
uqvect, vqvect = mpcalc.q_vector(u, v, t, 850*units.hPa, dx, dy)

# Compute the divergence of the Q-vectors calculated above
q_div = -2*mpcalc.divergence(uqvect, vqvect, dx=dx, dy=dy)


# Smooth the vertical motion so it is comparable to the QG omega forcing
w = mpcalc.smooth_n_point(ds.w, n_smooth, N_smooth_q)


### Make the plot ##############################################
print('Plotting...')


u = ds.u
v = ds.v
t = ds.t
z = ds.z
msl = ds.msl
pt = ds.pt

# Set the map projection (how the data will be displayed)
mapcrs = ccrs.PlateCarree()

# Set the data projection (GFS is lat/lon format)
datacrs = ccrs.PlateCarree()

# Start the figure and set an extent to only display a smaller graphics area
fig = plt.figure(1, figsize=(14, 12))
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([lons[0], lons[1], lats[1], lats[0]], ccrs.PlateCarree())

# Add map features to plot coastlines and state boundaries
ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray")
ax.add_feature(cfeature.COASTLINE.with_scale("50m"))
ax.add_feature(cfeature.BORDERS.with_scale("50m"), linestyle=":")

# Add gridlines with nicely spaced labels
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {"size": 10}
gl.ylabel_style = {"size": 10}
gl.xlocator = plt.FixedLocator(range(int(np.floor(lons[0]/10)*10),int( np.ceil(lons[1]/10)*10), 10))  # every 10 degrees longitude
gl.ylocator = plt.FixedLocator(range(int(np.ceil(lats[1]/10)*10), int(np.floor(lats[0]/10)*10), 10))   # every 10 degrees latitude


# Create colormap: blue to red for the divergence
cmap = plt.cm.bwr


# Plot 850-hPa Q-Vector Divergence and scale
clevs_qdiv = list(range(-20, -3, 4))+list(range(4, 21, 4))
cf = ax.contourf(ds.longitude, ds.latitude, q_div*1e18, clevs_qdiv, cmap=cmap,
                 extend='both', transform=datacrs)
cb = plt.colorbar(cf, orientation='horizontal', pad=0.05, aspect=50, extendrect=True,
                  ticks=clevs_qdiv)
cb.set_label('Q-Vector Div. (*10$^{18}$ m s$^{-1}$ kg$^{-1}$)')

# Plot 850-hPa Temperatures
clevs_t = np.arange(230, 310, 2)
csf = ax.contour(ds.longitude, ds.latitude, t, clevs_t, colors='grey',
                 linestyles='dashed', transform=datacrs)
plt.clabel(csf, fmt='%d')


# Plot 850-hPa Geopotential Heights
clevs_850_hght = np.arange(0, 8000, 30)
cs = ax.contour(ds.longitude, ds.latitude, z/9.81, clevs_850_hght, colors='black', transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot 850-hPa w
clev_w = [0.1, 0.15, 0.2]
cs = ax.contour(ds.longitude, ds.latitude, -w, clev_w, cmap='summer', transform=datacrs)


# Plot 850-hPa Q-vectors, scale to get nice sized arrows
wind_slice = (slice(None, None, 2), slice(None, None, 2))
ax.quiver(ds.longitude[wind_slice[0]], ds.latitude[wind_slice[1]],
          uqvect[wind_slice],
          vqvect[wind_slice],
          pivot='mid', color='black',
          scale=0.5e-11, scale_units='inches',
          transform=datacrs)

# Add some titles
plt.title('ERA5 850-hPa: z (black; m), T (grey; K)'
          ', Q-Vec. (arrows), Q-vec div. (colors), ascent (green)', loc='left',fontsize=10)

time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
plt.title('{} UTC'.format(time_str), loc='right',fontsize=10)

time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")

plt.savefig("Figures/ERA5_Q-vectors_"+time_str+".png", dpi=300, bbox_inches="tight")
#plt.show()
plt.close(fig)



### Plot the second figure #########################################

# Set the map projection (how the data will be displayed)
mapcrs = ccrs.PlateCarree()

# Set the data projection (GFS is lat/lon format)
datacrs = ccrs.PlateCarree()

# Start the figure and set an extent to only display a smaller graphics area
fig = plt.figure(1, figsize=(14, 12))
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([lons[0], lons[1], lats[1], lats[0]], ccrs.PlateCarree())

# Add map features to plot coastlines and state boundaries
ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray")
ax.add_feature(cfeature.COASTLINE.with_scale("50m"))
ax.add_feature(cfeature.BORDERS.with_scale("50m"), linestyle=":")

# Add gridlines with nicely spaced labels
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {"size": 10}
gl.ylabel_style = {"size": 10}
gl.xlocator = plt.FixedLocator(range(int(np.floor(lons[0]/10)*10),int( np.ceil(lons[1]/10)*10), 10))  # every 10 degrees longitude
gl.ylocator = plt.FixedLocator(range(int(np.ceil(lats[1]/10)*10), int(np.floor(lats[0]/10)*10), 10))   # every 10 degrees latitude


# Create colormap: blue to red for the divergence
cmap = plt.cm.nipy_spectral


# Plot tropopause potential temperature
clevs_pt = list(range(285, 385, 5))
cf = ax.contourf(ds_orig.longitude, ds_orig.latitude, ds_orig.pt, clevs_pt, cmap=cmap,
                 extend='both', transform=datacrs)
cbp = plt.colorbar(cf, orientation='horizontal', pad=0.05, aspect=50, extendrect=True,
                  ticks=clevs_pt)
cbp.set_label('Theta on DT. (K)')


# Plot Mean sea-level pressure
clevs_msl = np.arange(920, 1040, 4)
cs = ax.contour(ds.longitude, ds.latitude, msl/100, clevs_msl, colors='black', transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot 850-hPa Temperatures
clevs_t = np.arange(230, 310, 2)
csf = ax.contour(ds.longitude, ds.latitude, t, clevs_t, cmap='gist_rainbow',
                  transform=datacrs)
plt.clabel(csf, fmt='%d')


# Add some titles
plt.title('ERA5: MSLP (black; hPa), 850hPa T (K),theta on DT (colors; K)', loc='left',fontsize=10)

time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M")
plt.title('{} UTC'.format(time_str), loc='right',fontsize=10)

time_str = pd.Timestamp(ds.time.values).strftime("%Y-%m-%d_%H:%M")

plt.savefig("Figures/ERA5_tropopause_"+time_str+".png", dpi=300, bbox_inches="tight")
plt.show()



