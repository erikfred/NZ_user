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

g = 9.81
pad = 36 # this trims the ends after the low pass so there are no nan's

ncoutdir = '../LO_output/allinone/pickles/'

#### end setup ####

# LOADING
t_arr = pickle.load(open((ncoutdir + 't_arr.p'), 'rb'))
bp_tot = pickle.load(open((ncoutdir + 'bp_tot.p'), 'rb'))
bp_tot2 = pickle.load(open((ncoutdir + 'bp_tot2.p'), 'rb'))
bp_bc = pickle.load(open((ncoutdir + 'bp_bc.p'), 'rb'))
bp_ssh = pickle.load(open((ncoutdir + 'bp_ssh.p'), 'rb'))
lat = pickle.load(open((ncoutdir + 'lat.p'), 'rb'))
lon = pickle.load(open((ncoutdir + 'lon.p'), 'rb'))
bath = pickle.load(open((ncoutdir + 'bath.p'), 'rb'))

# convert to anomalies
bpm = np.mean(bp_tot, axis=0)
bp_anom = bp_tot - bpm; del bp_tot
bpm = np.mean(bp_tot2, axis=0)
bp_anom2 = bp_tot2 - bpm; del bp_tot2
bpm = np.mean(bp_bc, axis=0)
bp_bc_anom = bp_bc - bpm; del bp_bc
bpm = np.mean(bp_ssh, axis=0)
bp_ssh_anom = bp_ssh - bpm; del bp_ssh

ssh_anom = bp_ssh_anom / g / 1025

# print(ssh_anom.shape) # (48, 975, 270)

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.thermal

plt.close('all')

# strike-parallel line
isb = np.argwhere((bath<302) & (bath>298))
for nn in range(int(np.around(len(isb)/5))):
    # low pass filter
    anomlp = zfun.lowpass(bp_anom[:,isb[nn*5][0],isb[nn*5][1]], f='godin')[pad:-pad:24]
    anom2lp = zfun.lowpass(bp_anom2[:,isb[nn*5][0],isb[nn*5][1]], f='godin')[pad:-pad:24]
    bclp = zfun.lowpass(bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]], f='godin')[pad:-pad:24]
    sshlp = zfun.lowpass(bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]], f='godin')[pad:-pad:24]
    if nn==0:
        tlpd = t_arr[pad:-pad:24]
        # tlpd = []
        # for tt in range(len(tlp)):
        #     tlpd[tt] = datetime.fromtimestamp(tlp[tt])

    fig0 = plt.figure(figsize=(8,8))
    ax0 = fig0.add_subplot(111)
    # ax1.plot(tlp,anomlp,label='anom')
    # ax1.plot(tlpd,anom2lp,label='anom2')
    ax0.plot(tlpd,bclp,label='baroclinic')
    ax0.plot(tlpd,sshlp,label='ssh')
    ax0.plot(tlpd,bclp+sshlp,color='r',label='sum')

    if nn==0:
        yl1, yl2 = ax0.get_ylim()
    ax0.set_ylim(yl1-500,yl2+500)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1.5)
    ax0.grid(True)

    if not os.path.exists(ncoutdir + 'bp_plots'):
        os.mkdir(ncoutdir + 'bp_plots')

    # plt.show()
    plt.savefig(ncoutdir + 'bp_plots/mooring2_' + str(nn))
    plt.close()

"""
# strike-perpendicular line
for nn in range(20):
    # low pass filter
    anomlp = zfun.lowpass(bp_anom[:,30,nn*10], f='godin')[pad:-pad:24]
    anom2lp = zfun.lowpass(bp_anom2[:,30,nn*10], f='godin')[pad:-pad:24]
    bclp = zfun.lowpass(bp_bc_anom[:,30,nn*10], f='godin')[pad:-pad:24]
    sshlp = zfun.lowpass(bp_ssh_anom[:,30,nn*10], f='godin')[pad:-pad:24]
    if nn==0:
        tlpd = t_arr[pad:-pad:24]
        # tlpd = []
        # for tt in range(len(tlp)):
        #     tlpd[tt] = datetime.fromtimestamp(tlp[tt])

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)
    # ax1.plot(tlp,anomlp,label='anom')
    # ax1.plot(tlpd,anom2lp,label='anom2')
    ax1.plot(tlpd,bclp,label='baroclinic')
    ax1.plot(tlpd,sshlp,label='ssh')
    ax1.plot(tlpd,bclp+sshlp,color='r',label='sum')

    if nn==0:
        yl1, yl2 = ax1.get_ylim()
    ax1.set_ylim(yl1-500,yl2+500)
    ax1.set_title(str(np.around(lat[30,nn*10],decimals=2)) + ', ' + str(np.around(lon[30,nn*10],decimals=2))
        + ', ' + str(np.around(bath[30,nn*10])))
    ax1.legend()
    ax1.axhline(c='k',lw=1.5)
    ax1.grid(True)

    if not os.path.exists(ncoutdir + 'bp_plots'):
        os.mkdir(ncoutdir + 'bp_plots')

    # plt.show()
    plt.savefig(ncoutdir + 'bp_plots/mooring_' + str(nn))
    plt.close()
"""

# # map view plots
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
