"""
plot_components.py
Generate time series plots of baroclinic, sea surface height, and total
bottom pressure contributions.
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import cmocean
import matplotlib.pyplot as plt
import pickle

from lo_tools import Lfun, zfun, zrfun

g = 9.81
savedir = '../LO_output/timeseries/bp_components/'
topdir = '../LO_output/allinone/'
loadir = topdir + 'pickles_2017-18/'

#### end setup ####

# LOADING
tlp = pickle.load(open((loadir + 'tlp.p'), 'rb'))
bp_tot = pickle.load(open((loadir + 'bp_tot.p'), 'rb'))
bp_tot2 = pickle.load(open((loadir + 'bp_tot2.p'), 'rb'))
bp_bc = pickle.load(open((loadir + 'bp_bc.p'), 'rb'))
bp_ssh = pickle.load(open((loadir + 'bp_ssh.p'), 'rb'))
lat = pickle.load(open((loadir + 'lat.p'), 'rb'))
lon = pickle.load(open((loadir + 'lon.p'), 'rb'))
bath = pickle.load(open((loadir + 'bath.p'), 'rb'))

# CONVERT TO ANOMALIES
# annual means
Tot = np.mean(bp_tot, axis=0)
Tot2 = np.mean(bp_tot2, axis=0)
Bc = np.mean(bp_bc, axis=0)
Ssh = np.mean(bp_ssh, axis=0)
# anomalies
bp_anom = bp_tot - Tot
bp_anom2 = bp_tot2 - Tot2
bp_bc_anom = bp_bc - Bc
bp_ssh_anom = bp_ssh - Ssh

ssh_anom = bp_ssh_anom / g / 1025

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
    fig0 = plt.figure(figsize=(8,8))
    ax0 = fig0.add_subplot(311)
    # ax1.plot(tlp,anomlp,label='anom')
    # ax1.plot(tlp,anom2lp,label='anom2')
    ax0.plot(tlp,bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]],label='baroclinic')
    ax0.plot(tlp,bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]],label='ssh')
    ax0.plot(tlp,bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]+bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]],color='r',label='sum')

    if nn==0:
        yl1, yl2 = ax0.get_ylim()
    ax0.set_ylim(yl1-500,yl2+500)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1.5)
    ax0.grid(True)

    ax1 = fig0.add_subplot(312)
    x1 = np.convolve(bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]], np.ones(10)/10, mode='same')
    x2 = np.convolve(bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]], np.ones(10)/10, mode='same')
    metric = x1/x2
    ax1.plot(tlp,metric,label='bc/ssh')
    ax1.legend()
    ax1.set_ylim(-2,2)
    ax1.axhline(c='k',lw=1.5)
    ax1.grid(True)

    ax2 = fig0.add_subplot(313)
    ax2.plot(tlp,bp_anom[:,isb[nn*5][0],isb[nn*5][1]],label='anom')
    ax2.plot(tlp,bp_anom2[:,isb[nn*5][0],isb[nn*5][1]],label='anom2')
    ax2.legend()
    ax2.set_ylim(yl1-500,yl2+500)
    ax2.axhline(c='k',lw=1.5)
    ax2.grid(True)

    if not os.path.exists(savedir + 'strike_par'):
        os.mkdir(savedir + 'strike_par')

    # plt.show()
    plt.savefig(savedir + 'strike_par/mooring_' + str(nn))
    plt.close()

# strike-perpendicular line
for nn in range(20):
    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(311)
    ax1.plot(tlp,bp_bc_anom[:,30,nn*10],label='baroclinic')
    ax1.plot(tlp,bp_ssh_anom[:,30,nn*10],label='ssh')
    ax1.plot(tlp,bp_bc_anom[:,30,nn*10]+bp_ssh_anom[:,30,nn*10],color='r',label='sum')

    if nn==0:
        yl1, yl2 = ax1.get_ylim()
    ax1.set_ylim(yl1-500,yl2+500)
    ax1.set_title(str(np.around(lat[30,nn*10],decimals=2)) + ', ' + str(np.around(lon[30,nn*10],decimals=2))
        + ', ' + str(np.around(bath[30,nn*10])))
    ax1.legend()
    ax1.axhline(c='k',lw=1.5)
    ax1.grid(True)

    ax2 = fig1.add_subplot(312)
    x1 = np.convolve(bp_bc_anom[:,30,nn*10], np.ones(10)/10, mode='same')
    x2 = np.convolve(bp_ssh_anom[:,30,nn*10], np.ones(10)/10, mode='same')
    metric = x1/x2
    ax2.plot(tlp,metric,label='bc/ssh')
    ax2.legend()
    ax2.set_ylim(-2,2)
    ax2.axhline(c='k',lw=1.5)
    ax2.grid(True)

    ax3 = fig1.add_subplot(313)
    ax3.plot(tlp,bp_anom[:,30,nn*10],label='anom')
    ax3.plot(tlp,bp_anom2[:,30,nn*10],label='anom2')
    ax3.legend()
    ax3.set_ylim(yl1-500,yl2+500)
    ax3.axhline(c='k',lw=1.5)
    ax3.grid(True)

    if not os.path.exists(savedir + 'strike_perp'):
        os.mkdir(savedir + 'strike_perp')

    # plt.show()
    plt.savefig(savedir + 'strike_perp/mooring_' + str(nn))
    plt.close()
