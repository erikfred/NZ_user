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
topdir = '../LO_output/allinone/'
loadir = topdir + 'pickles_2017-18/'
outdir = '../LO_output/timeseries/'
savedir = outdir + loadir[-8:]

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

tlp2 = [datetime.fromtimestamp(t) for t in tlp]

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 1 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.balance

plt.close('all')

# PRESSURE COMPONENTS
# strike-parallel line at 200 m
isb = np.argwhere((bath<202) & (bath>198))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/10))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(323)
    ax0.plot(tlp2,bp_bc_anom[:,isb[nn*10][0],isb[nn*10][1]]/100,lw=lw,label='steric')
    ax0.plot(tlp2,bp_ssh_anom[:,isb[nn*10][0],isb[nn*10][1]]/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,bp_anom2[:,isb[nn*10][0],isb[nn*10][1]]/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*10][0],isb[nn*10][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*10][0],isb[nn*10][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*10][0],isb[nn*10][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[10][0],isb[10][1]],lon[isb[10][0],isb[10][1]],'^r',markersize=mk-5)
    ax1.plot(lat[isb[nn*10][0],isb[nn*10][1]],lon[isb[nn*10][0],isb[nn*10][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_200'):
        os.mkdir(savedir + 'strike_par_200')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_200/mooring_' + str(nn))
    plt.savefig(savedir + 'strike_par_200/mooring_' + str(nn) + '.eps')
    plt.close()

# strike-parallel line at 300 m
isb = np.argwhere((bath<302) & (bath>298))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/5))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(323)
    ax0.plot(tlp2,bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,label='steric')
    ax0.plot(tlp2,bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[stalist,0],isb[stalist,1]],lon[isb[stalist,0],isb[stalist,1]],'ok',markersize=mk-5)
    ax1.plot(lat[isb[nn*5][0],isb[nn*5][1]],lon[isb[nn*5][0],isb[nn*5][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_300'):
        os.mkdir(savedir + 'strike_par_300')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_300/mooring_' + str(nn))
    plt.savefig(savedir + 'strike_par_300/mooring_' + str(nn) + '.eps')
    plt.close()

# strike-parallel line at 2000 m
isb = np.argwhere((bath<2002) & (bath>1998))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/5))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(323)
    ax0.plot(tlp2,bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,label='steric')
    ax0.plot(tlp2,bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[stalist,0],isb[stalist,1]],lon[isb[stalist,0],isb[stalist,1]],'ok',markersize=mk-5)
    ax1.plot(lat[isb[nn*5][0],isb[nn*5][1]],lon[isb[nn*5][0],isb[nn*5][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_2000'):
        os.mkdir(savedir + 'strike_par_2000')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_2000/mooring_' + str(nn))
    plt.savefig(savedir + 'strike_par_2000/mooring_' + str(nn) + '.eps')
    plt.close()

# # strike-perpendicular line
# stalist = [x * 10 for x in list(range(5,24))]
# for nn in range(5,24):
#     for mm in range(1,5):
#         fig1 = plt.figure(figsize=(8,8))
#         ax1 = fig1.add_subplot(212)
#         ax1.plot(tlp2,bp_bc_anom[:,150*mm-50,nn*10],label='steric')
#         ax1.plot(tlp2,bp_ssh_anom[:,150*mm-50,nn*10],label='eustatic')
#         ax1.plot(tlp2,bp_bc_anom[:,150*mm-50,nn*10]+bp_ssh_anom[:,150*mm-50,nn*10],color='r',label='total')
#
#         ax1.set_ylim(-1500,1500)
#         ax1.set_title(str(np.around(lat[150*mm-50,nn*10],decimals=2)) + ', ' + str(np.around(lon[150*mm-50,nn*10],decimals=2))
#             + ', ' + str(np.around(bath[150*mm-50,nn*10])))
#         ax1.legend()
#         ax1.axhline(c='k',lw=1)
#         ax1.grid(True)
#
#         ax2 = fig1.add_subplot(321)
#         bth = ax2.contour(lat, lon, bath, [4, 300, 2000], colors='black')
#         ax2.plot(lat[150*mm-50,stalist],lon[150*mm-50,stalist],'ok',markersize=mk-5)
#         ax2.plot(lat[150*mm-50,nn*10],lon[150*mm-50,nn*10],'or',markersize=mk-5)
#         ax2.invert_yaxis()
#         ax2.grid(True)
#
#         if not os.path.exists(savedir + 'strike_perp_' + str(mm)):
#             os.mkdir(savedir + 'strike_perp_' + str(mm))
#
#         # plt.show()
#         plt.savefig(savedir + 'strike_perp_' + str(mm) + '/mooring_' + str(nn-5))
#         plt.close()

# DIFFERENCES RELATIVE TO REFERENCE LOCATION
# strike-parallel line at 200 m
isb = np.argwhere((bath<202) & (bath>198))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/10))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(324)
    ax0.plot(tlp2,(bp_bc_anom[:,isb[nn*10][0],isb[nn*10][1]]-bp_bc_anom[:,isb[10][0],isb[10][1]])/100,lw=lw,label='steric')
    ax0.plot(tlp2,(bp_ssh_anom[:,isb[nn*10][0],isb[nn*10][1]]-bp_ssh_anom[:,isb[10][0],isb[10][1]])/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,(bp_anom2[:,isb[nn*10][0],isb[nn*10][1]]-bp_anom2[:,isb[10][0],isb[10][1]])/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*10][0],isb[nn*10][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*10][0],isb[nn*10][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*10][0],isb[nn*10][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[10][0],isb[10][1]],lon[isb[10][0],isb[10][1]],'^r',markersize=mk-5)
    ax1.plot(lat[isb[nn*10][0],isb[nn*10][1]],lon[isb[nn*10][0],isb[nn*10][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_200'):
        os.mkdir(savedir + 'strike_par_200')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_200/difference_' + str(nn))
    plt.savefig(savedir + 'strike_par_200/difference_' + str(nn) + '.eps')
    plt.close()

# strike-parallel line at 300 m
isb = np.argwhere((bath<302) & (bath>298))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/5))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(324)
    ax0.plot(tlp2,(bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]-bp_bc_anom[:,isb[0][0],isb[0][1]])/100,lw=lw,label='steric')
    ax0.plot(tlp2,(bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]-bp_ssh_anom[:,isb[0][0],isb[0][1]])/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,(bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]-bp_anom2[:,isb[0][0],isb[0][1]])/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[0][0],isb[0][1]],lon[isb[0][0],isb[0][1]],'^r',markersize=mk-5)
    ax1.plot(lat[isb[nn*5][0],isb[nn*5][1]],lon[isb[nn*5][0],isb[nn*5][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_300'):
        os.mkdir(savedir + 'strike_par_300')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_300/difference_' + str(nn))
    plt.savefig(savedir + 'strike_par_300/difference_' + str(nn) + '.eps')
    plt.close()

# strike-parallel line at 2000 m
isb = np.argwhere((bath<2002) & (bath>1998))
stalist = [x*5 for x in list(range(int(np.around(len(isb)/5))))]
for nn in range(int(np.around(len(isb)/5))):
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(324)
    ax0.plot(tlp2,(bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]-bp_bc_anom[:,isb[0][0],isb[0][1]])/100,lw=lw,label='steric')
    ax0.plot(tlp2,(bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]-bp_ssh_anom[:,isb[0][0],isb[0][1]])/100,lw=lw,label='eustatic')
    ax0.plot(tlp2,(bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]-bp_anom2[:,isb[0][0],isb[0][1]])/100,lw=lw,color='r',label='total')

    ax0.set_ylim(-15,15)
    ax0.set_title(str(np.around(lat[isb[nn*5][0],isb[nn*5][1]],decimals=2)) +
        ', ' + str(np.around(lon[isb[nn*5][0],isb[nn*5][1]],decimals=2))
        + ', ' + str(np.around(bath[isb[nn*5][0],isb[nn*5][1]])))
    ax0.legend()
    ax0.axhline(c='k',lw=1)
    ax0.grid(True)
    ax0.tick_params(labelsize=fs-2)
    plt.xticks(rotation=45)

    ax1 = fig0.add_subplot(321)
    bth = ax1.contour(lat, lon, bath, [4, 300, 2000], colors='black')
    ax1.plot(lat[isb[0][0],isb[0][1]],lon[isb[0][0],isb[0][1]],'^r',markersize=mk-5)
    ax1.plot(lat[isb[nn*5][0],isb[nn*5][1]],lon[isb[nn*5][0],isb[nn*5][1]],'or',markersize=mk-5)
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + 'strike_par_2000'):
        os.mkdir(savedir + 'strike_par_2000')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + 'strike_par_2000/difference_' + str(nn))
    plt.savefig(savedir + 'strike_par_2000/difference_' + str(nn) + '.eps')
    plt.close()
