import xarray as xr
import numpy as np
import gsw
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

# set mooring extraction to analyze (just needs salt, temp, and zeta)
Ldir = Lfun.Lstart()

out_dir = Ldir['parent'] / 'LO_output' / 'bpress_PM' # I think this assumes I'm working out of the bpress subfolder
Lfun.make_dir(out_dir)

# sn_list = ['CE01', 'CE02', 'CE04', 'PN01A'] # OOI cabled stations?
sn_list = ['PN01A']

fs = 12 # fontsize
def icb(ax, cs):
    # cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2)
    # cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    ax.fill([.93,.995,.995,.93],[.05,.05,.95,.95],'w', alpha=.5, transform=ax.transAxes)
    cbaxes = inset_axes(ax, width="2%", height="80%", loc='right', borderpad=3)
    cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    cb.ax.tick_params(labelsize=.85*fs)

plt.close('all')
for sn in sn_list:

    fn = Ldir['LOo'] / 'extract' / 'cas6_v3_lo8b' / 'moor' / 'ooi' / (sn + '_2018.01.01_2018.12.31.nc')
    # ‘/data1/parker/LO_roms/cas6_v0_live’
    print(fn)
