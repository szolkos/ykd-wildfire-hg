import numpy as np
import pandas as pd
import xarray as xr
import cmocean
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib as mpl

import json
with open('../misc/sargassum.json', 'r') as f:
    sargassum = LinearSegmentedColormap("sargassum", json.load(f))

with open('../misc/nature_reviews.json', 'r') as f:
    nature = json.load(f)

def get_domain_bounds():
    return {'Alaska':{'lat_min': 51,    'lat_max': 73,   'lon_min': -170,     'lon_max': -138.5,},
            'YKD':   {'lat_min': 58,    'lat_max': 64,   'lon_min': -168,     'lon_max': -160,},
            'IKU':   {'lat_min': 61.25, 'lat_max': 61.5, 'lon_min': -163.125, 'lon_max': -162.5,},
            }

domain_bounds = get_domain_bounds()
# ------------------------------------------------------------------------
# Set up general plotting parameters
IKU_fire_lat_min, IKU_fire_lat_max = domain_bounds['IKU']['lat_min'], domain_bounds['IKU']['lat_max']
IKU_fire_lon_min, IKU_fire_lon_max = domain_bounds['IKU']['lon_min'], domain_bounds['IKU']['lon_max']

YKD_lat_min, YKD_lat_max = domain_bounds['YKD']['lat_min'], domain_bounds['YKD']['lat_max']
YKD_lon_min, YKD_lon_max = domain_bounds['YKD']['lon_min'], domain_bounds['YKD']['lon_max']

Domain_lat_min, Domain_lat_max = domain_bounds['Alaska']['lat_min'], domain_bounds['Alaska']['lat_max']
Domain_lon_min, Domain_lon_max = domain_bounds['Alaska']['lon_min'], domain_bounds['Alaska']['lon_max']

Plot_domain_lat_min, Plot_domain_lat_max = Domain_lat_min+1, Domain_lat_max-1
Plot_domain_lon_min, Plot_domain_lon_max = Domain_lon_min+1, Domain_lon_max-1
# ------------------------------------------------------------------------

def add_bounding_box(ax, lat_min, lat_max, lon_min, lon_max, lat_res=0.25, lon_res=0.3125, **kwargs):
    ''' plot box around region of interest. lat and lon are the center of the box.
        lat_res and lon_res are the resolution of the box.
        Function will plot a box with the specified resolution around the cells contained by the
        grid centers.'''

    ax.add_patch(mpatches.Rectangle(xy=[(lon_min-0.5*lon_res), (lat_min-0.5*lat_res)], 
                                    width=(lon_max-lon_min+lon_res), height=(lat_max-lat_min+lat_res),
                                    facecolor='none', transform=ccrs.PlateCarree(), **kwargs))
    return

def convert_ppqv_to_ngm3(ds):
    ''' convert ppqv to ng m-3 '''
    return ds*200.59*44.642*1e9

def plot_data2(ax, data, title='', cmap=sargassum, levels=None, vmin=None, vmax=None, norm=None,
               lon_min=Plot_domain_lon_min, lon_max=Plot_domain_lon_max, lat_min=Plot_domain_lat_min, lat_max=Plot_domain_lat_max,
               contour=False, y_gridlines=[55, 60, 65, 70], x_gridlines=[-165, -155, -145], coastline_lw=0.5, gridline_lw=0.5):

    # The projection of the map
    noProj = ccrs.PlateCarree()

    # This plots parallel and meridian arcs around a target area that will be used as the map boundary
    [ax_hdl] = ax.plot([lon_min, lon_max, lon_max, lon_min, lon_min],
                       [lat_min, lat_min, lat_max, lat_max, lat_min],
                        color='black', linewidth=0.5, marker='none', transform=noProj, zorder=5)
    # Get the `Path` of the plot
    tx_path = ax_hdl._get_transformed_path()
    path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()

    # Use the path's vertices to create a polygon
    polygon = mpath.Path( path_in_data_coords.vertices)
    ax.set_boundary(polygon) #This masks-out unwanted part of the plot
    
    ax.add_feature(cartopy.feature.COASTLINE, edgecolor='black', lw=coastline_lw, zorder=4)
    
    # -- plot the data
    if contour:
        im = data.plot.contourf(ax=ax, transform=noProj, vmin=vmin, vmax=vmax,
                               cmap=cmap, levels=levels, add_colorbar=False, zorder=2)
    else:
        im = data.plot(ax=ax, transform=noProj, vmin=vmin, vmax=vmax, norm=norm,
                       cmap=cmap, levels=levels, add_colorbar=False, zorder=2)
        
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False, linewidth=gridline_lw, alpha=0.2, linestyle='-')

    # -- add ticks
    gl.ylocator = mticker.FixedLocator(y_gridlines)
    gl.xlocator = mticker.FixedLocator(x_gridlines)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': '0.3'}
    gl.ylabel_style = {'size': 8, 'color': '0.3'}

    add_bounding_box(ax, YKD_lat_min, YKD_lat_max, YKD_lon_min, YKD_lon_max, edgecolor='k', lw=0.5, zorder=5)
    ax.set_title(title, pad=0.1)

    return ax, im

def plot_inset(ds, ax, contour=False, ocean_color='None', plot_rivers=False, **kwargs):

    if contour:
        im = ds.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), **kwargs)
    else:
        im = ds.plot(ax=ax, transform=ccrs.PlateCarree(), **kwargs)

    ax.add_feature(cfeature.BORDERS, linestyle='-', lw=1, color='k', alpha=0.2,)
    ax.add_feature(cfeature.COASTLINE, linestyle='-', lw=1, color='k', alpha=0.5,)
    ax.add_feature(cartopy.feature.OCEAN, edgecolor='None', facecolor=ocean_color, zorder=2)
    if plot_rivers:
        ax.add_feature(cartopy.feature.RIVERS, zorder=1, lw=0.5, edgecolor='k', alpha=0.2)
    
    ax.set_extent([YKD_lon_min, YKD_lon_max, YKD_lat_min, YKD_lat_max], crs=ccrs.PlateCarree())
    # add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, alpha=0.2, linestyle='-',)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': '0.3'}
    gl.ylabel_style = {'size': 8, 'color': '0.3'}

    return ax, im


# colormap utilities

# convert 0 - 255 rgb to 0 - 1
def rgb_to_1(rgb):
    return [c / 255 for c in rgb]

# create matplotlib colormap from rgb
def create_colormap(name, colors):
    colors = [rgb_to_1(c) for c in colors]
    cm = mpl.colors.ListedColormap([tuple(c) for c in colors])
    cm.name = name
    return cm

# create matplotlib colormap from json file
def create_colormap_from_json(json_file):
    with open(json_file) as f:
        data = json.load(f)

    name = data['name']
    colors = data['colors']
    return create_colormap(name, colors)

# create matplotlib colormap from json file
def create_colormap_from_json_file(json_file):
    with open(json_file) as f:
        data = json.load(f)

    name = data['name']
    colors = data['colors']
    return create_colormap(name, colors)

# display all colormaps
def display_colormaps(colormaps):
    n = len(colormaps)
    fig, axs = plt.subplots(n, 1, figsize=(6, n * 2))
    for i, cm in enumerate(colormaps):
        plt.sca(axs[i])
        # create an n, 1 array to display the colormap
        x = np.linspace(0, 1, len(cm.colors)).reshape(1, -1)
        plt.imshow(x, cmap=cm, aspect=0.7)
        plt.xticks([])
        plt.yticks([])
        plt.title(cm.name)
    # remove vertical space between subplots
    fig.subplots_adjust(hspace=-0.8)
    plt.show()