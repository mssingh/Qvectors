# map_utils.py
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import geopandas as gpd


melb_shape_file = './SA4_2021_AUST_GDA94/SA4_2021_AUST_GDA94.shp'


def make_map(lons,lats,grid_spacing=10,states=False,Melbourne=False):


    ## Figure ###########################################################################

    # Set the map projection
    mapcrs = ccrs.PlateCarree()

    # Start the figure and set an extent to only display a smaller graphics area
    fig = plt.figure(figsize=(17, 12))
    ax = plt.subplot(111, projection=mapcrs)
    ax.set_extent([lons[0], lons[1], lats[1], lats[0]], ccrs.PlateCarree())


    ## Boundaries #######################################################################

    # Coastlines
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),color='black',linewidth=3)


    if states:

        # State boundaries
        states_provinces = cfeature.NaturalEarthFeature(
           category='cultural',
           name='admin_1_states_provinces_lines',
           scale='50m',
           facecolor='none'
        )
        ax.add_feature(states_provinces, edgecolor='black', linewidth=1)


    if Melbourne:

        # Plot the Greater Melbourne boundary
        gccsa = gpd.read_file(melb_shape_file)
        melb_metro = gccsa[gccsa['GCC_NAME21'] == 'Greater Melbourne']
        melb_outline = melb_metro.dissolve()

        ax.add_geometries(melb_outline.geometry, crs=ccrs.PlateCarree(),
              facecolor='none', edgecolor='red', linewidth=1)

        # PLot the location of Monash
        ax.plot(145.1339,-37.9078, marker='o', color='red', markersize=8, transform=ccrs.PlateCarree())


    ## Gridlines

    G = grid_spacing

    # Add gridlines with nicely spaced labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
    gl.top_labels = False
    gl.right_labels = True
    gl.xlabel_style = {"size": 14}
    gl.ylabel_style = {"size": 14}
    gl.xlocator = plt.FixedLocator(range(int(np.floor(lons[0]/G)*G),int( np.ceil(lons[1]/G)*G), G))  # every G degrees longitude
    gl.ylocator = plt.FixedLocator(range(int(np.floor(lats[0]/G)*G), int(np.ceil(lats[1]/G)*G), G))   # every G degrees latitude



    return fig,ax
