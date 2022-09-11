"""
bp_plotter.py
Loads existing bp data from pickle files and makes exploratory plots
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import cmocean
import matplotlib.pyplot as plt
import pickle

sys.path.append(os.path.abspath('util'))
import Lfun
import zrfun
import zfun

ncoutdir = '../LO_output/spinup/bp_extractions/'

#### end setup ####

# LOADING
t_arr = pickle.load(open((ncoutdir + 't_arr.p'), 'rb'))
bp_anom = pickle.load(open((ncoutdir + 'bp_anom.p'), 'rb'))
bp_anom2 = pickle.load(open((ncoutdir + 'bp_anom2.p'), 'rb'))
bp_bc_anom = pickle.load(open((ncoutdir + 'bp_bc.p'), 'rb'))
bp_ssh_anom = pickle.load(open((ncoutdir + 'bp_ssh.p'), 'rb'))
lat = pickle.load(open((ncoutdir + 'lat.p'), 'rb'))
lon = pickle.load(open((ncoutdir + 'lon.p'), 'rb'))

# print(ssh_anom.shape) # (48, 975, 270)

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.thermal

plt.close('all')

fig1 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(111)
# pts = ax1.plot(t_arr[0:58],bp_anom[0:58,30,100],label='anom')
# ax1.plot(t_arr[0:58],bp_anom2[0:58,30,100],label='anom2')
ax1.plot(t_arr[0:58],bp_bc_anom[0:58,30,100],label='baroclinic')
# ax1.plot(t_arr[0:58],bp_ssh_anom[0:58,30,100],label='ssh')

ax1.set_title(str(lat[30,100]) + ', ' + str(lon[30,100]))
ax1.legend()

if not os.path.exists(ncoutdir + 'bp_maps'):
    os.mkdir(ncoutdir + 'bp_maps')

plt.show()
# plt.savefig(ncoutdir + 'bp_maps/mooring')
plt.close()

# for sv in range(nf): # for each day
#     fig = plt.figure(figsize=(6,10))
#     ax = fig.add_subplot(111)
#     cs = ax.pcolormesh(lon, lat, ssh_anom[sv,:,:], cmap=cmap, vmin=-0.2, vmax=0.2)
#     bth = ax.contour(lon, lat, bath, [300, 2000], colors='black')
#
#     ax.axis('square')
#     # ax.set_xlim((minlon,maxlon))
#     # ax.set_ylim((minlat, maxlat))
#     ax.grid(True)
#
#     cb = fig.colorbar(cs)
#     cb.set_label('SSH (m)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(ncoutdir + 'SSH_maps'):
#         os.mkdir(ncoutdir + 'SSH_maps')
#     plt.savefig(ncoutdir + 'SSH_maps/temp' + str(sv).zfill(4))
#     plt.close()
