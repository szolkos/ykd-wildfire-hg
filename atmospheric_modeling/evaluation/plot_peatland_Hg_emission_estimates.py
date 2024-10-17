# -- standalone script to plot peatland Hg emission estimates figure
import json
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from plotting_functions import *
import matplotlib.colors as mcolors
plt.style.use('../misc/acgc.mplstyle')

from Additional_Colormaps import *

if __name__ == "__main__":

    # load nature from json
    with open('../misc/nature_reviews.json') as f:
        nature = json.load(f)

    new_cmaps = load_all_cmaps()
    domains = get_domain_bounds()

    stones  = create_colormap('stone',  nature['stone'])
    greys   = create_colormap('grey',   nature['grey'])
    blues   = create_colormap('blue',   nature['blue'])

    # set axes.prop_cycle to custom color cycle
    new_colors = [stones.colors[4], greys.colors[5]]
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=new_colors)

    df = pd.read_csv('../misc/peat_emissions.csv')

    # - set ticklabel size
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    # set defult marker size
    plt.rcParams['lines.markersize'] = 4.5
    # set default marker edge width
    plt.rcParams['lines.markeredgewidth'] = 0.8

    c1 = '0.2'
    c2 = tuple([i/255 for i in nature['blue'][4]])

    v1, v1_lo, v1_hi = 'hg_em_med_peat', 'hg_em_lo_peat', 'hg_em_hi_peat'
    v2, v2_lo, v2_hi = 'hg_em_fxn_med_peat', 'hg_em_fxn_lo_peat', 'hg_em_fxn_hi_peat'

    # -- convert peatland fraction to percentage
    df[v2]    = df[v2]*100
    df[v2_lo] = df[v2_lo]*100
    df[v2_hi] = df[v2_hi]*100

    # -- set up figure
    fig = plt.figure(figsize=(3.25,2))
    ax = plt.subplot(111)

    # -- plot peatland wildfire Hg emissions in high-latitude NH
    ax.plot(df['year'], df[v1], marker='o', markeredgecolor=c1, markerfacecolor='w', c=c1, ls='-', zorder=2)
    ax.fill_between(df['year'], df[v1_lo], df[v1_hi], color=c1, alpha=0.1, lw=0, edgecolor='None', zorder=0)

    # -- plot peatland fraction of total high-latitude NH wildfire emissions
    ax2 = ax.twinx()
    ax2.plot(df['year'], df[v2], marker='o', markeredgecolor=c2, markerfacecolor='w', c=c2, ls='--', zorder=2)
    ax2.fill_between(df['year'], df[v2_lo], df[v2_hi], color=c2, alpha=0.1, lw=0, edgecolor='None', zorder=0)

    # -- set common y-axis params
    ymin, ymax = 0, 85
    yticks_major = [10, 20, 30, 40, 50, 60, 70, 80]
    yticks_minor = [5, 15, 25, 35, 45, 55, 65, 75]

    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=yticks_major)
    ax.set_yticks(ticks=yticks_minor, minor=True)

    ax2.set_ylim(ymin, ymax)
    ax2.set_yticks(ticks=yticks_major)
    ax2.set_yticks(ticks=yticks_minor, minor=True)

    # -- set x-axis params
    ax.set_xticks(ticks=[2005, 2010, 2015, 2020])
    ax.set_xticks(ticks=np.arange(2002, 2023, 1), minor=True)
    ax.set_xlim(2001, 2023)

    # change ax2 y-ticks to blue
    ax2.tick_params(axis='y', colors=c2)
    ax2.tick_params(axis='y', which='minor', colors=c2)

    ax.set_ylabel('Peatland Hg emissions (Mg y$^{-1}$)', fontsize=9)
    ax2.set_ylabel('Peatland fraction of total northern\nwildfire Hg emissions (%)', fontsize=9, rotation=270, labelpad=24, color=c2)

    ax.grid(False)
    ax2.grid(False)

    plt.tight_layout(pad=0)
    plt.savefig('../figures/misc/peatland_emissions.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('../figures/misc/peatland_emissions.png', format='png', bbox_inches='tight', dpi=2000)
    plt.close()