"""
pick_eddies.py
Generates daily plots of SSH from which user can select eddy centers to be saved
as a list of coordinates.
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import cmocean
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pickle
from time import mktime

from lo_tools import Lfun, zfun, zrfun

g = 9.81
topdir = '../LO_output/allinone/'
loadir = topdir + 'pickles_2017-18/'
outdir = '../LO_output/eddy_tracking/'
savedir = outdir + loadir[-8:]

# name the eddy
estr = 'eddy3';
pos = False; # set to True for + SSH anomaly and False for -
# save time by reducing this interval as much as possible
n1 = 195;
n2 = 227;

#### end setup ####

# plotting parameters
fs = 14 # primary fontsize
lw = 1 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.balance # formerly thermal
plt.close('all')

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

if os.path.isfile(savedir + estr + '/bp_eddy.p'):
    t_eddy = pickle.load(open((savedir + estr + '/t_eddy.p'), 'rb'))
    bp_anom_eddy = pickle.load(open((savedir + estr + '/bp_eddy.p'), 'rb'))
    bp_bc_eddy = pickle.load(open((savedir + estr + '/bc_eddy.p'), 'rb'))
    bp_ssh_eddy = pickle.load(open((savedir + estr + '/ssh_eddy.p'), 'rb'))
    ctrs_ll = pickle.load(open((savedir + estr + '/ctrs_ll.p'), 'rb'))
    ctrs_xy = pickle.load(open((savedir + estr + '/ctrs_xy.p'), 'rb'))

    t1 = mktime(t_eddy[0].timetuple())
    t2 = mktime(t_eddy[-1].timetuple())
    n1 = np.argmin(np.abs(tlp - t1))
    n2 = np.argmin(np.abs(tlp - t2)) + 1
else:
    # Track eddy and extract BP components at those points
    for nn in range(n2-n1): # for each day
        dn = n1 + nn
        di = datetime.fromtimestamp(tlp[dn])
        # ssh component
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        cs = ax.pcolormesh(lon, lat, bp_ssh_anom[dn,:,:], cmap=cmap, vmin=-1000, vmax=1000)
        bth = ax.contour(lon, lat, bath, [300, 2000], colors='black')

        ax.axis('square')
        ax.grid(True)
        ax.set_title('P_SSH ' + di.strftime('%m/%d/%y'))

        cb = fig.colorbar(cs)
        cb.set_label('Pressure (Pa)', fontsize=fs)
        cb.ax.tick_params(labelsize=fs)

        # GUI pick approximate center of eddy
        ctr=[];
        if nn==1:
            print('Select eddy center with mouse')
        ctr = np.asarray(plt.ginput(1, timeout=-1))
        # print(ctr)

        # use local extrema as true center
        ilo = np.argmin(np.abs(lon[0,:] - ctr[0,0]));
        ila = np.argmin(np.abs(lat[:,0] - ctr[0,1]));
        ssh_temp=bp_ssh_anom[dn,ila-15:ila+15,ilo-15:ilo+15].squeeze()
        if pos:
            ctr_alt = np.argwhere(ssh_temp == np.max(ssh_temp))
        else:
            ctr_alt = np.argwhere(ssh_temp == np.min(ssh_temp))
        # print(ctr_alt)

        ctr_alt2 = np.asarray([[ctr_alt[0,1], ctr_alt[0,0]]]);
        ctr_alt2[0,0] = ctr_alt2[0,0] + ilo - 15;
        ctr_alt2[0,1] = ctr_alt2[0,1] + ila - 15;
        # print(ctr_alt2)
        # print(lon[0,ctr_alt2[0,0]],lat[ctr_alt2[0,1],0])

        # # plot showing local SSH field and pick vs. relocation, for spot checking
        # fig = plt.figure(figsize=(8,8))
        # ax = fig.add_subplot(111)
        # ax.pcolormesh(lon[ila-15:ila+15,ilo-15:ilo+15],lat[ila-15:ila+15,ilo-15:ilo+15],ssh_temp, cmap=cmap)
        # ax.plot(ctr[0,0],ctr[0,1],'ko')
        # ax.plot(lon[0,ctr_alt2[0,0]],lat[ctr_alt2[0,1],0],'kx')
        # plt.show()

        if nn==0:
            ctrs_ll = np.asarray([[lon[0,ctr_alt2[0,0]],lat[ctr_alt2[0,1],0]]]);
            # print(ctrs_ll)
            ctrs_xy = ctr_alt2;
            # print(ctrs_xy)
            bp_anom_eddy = np.asarray([bp_anom2[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]);
            bp_bc_eddy = np.asarray([bp_bc_anom[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]);
            bp_ssh_eddy = np.asarray([bp_ssh_anom[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]);
            # print(bp_anom_eddy)
        else:
            ctrs_ll = np.append(ctrs_ll, np.asarray([[lon[0,ctr_alt2[0,0]],lat[ctr_alt2[0,1],0]]]), axis=0);
            # print(ctrs_ll)
            ctrs_xy = np.append(ctrs_xy, ctr_alt2, axis=0);
            # print(ctrs_xy)
            bp_anom_eddy = np.append(bp_anom_eddy, np.asarray([bp_anom2[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]), axis=0);
            bp_bc_eddy = np.append(bp_bc_eddy, np.asarray([bp_bc_anom[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]), axis=0);
            bp_ssh_eddy = np.append(bp_ssh_eddy, np.asarray([bp_ssh_anom[dn,ctr_alt2[0,1],ctr_alt2[0,0]]]), axis=0);

        plt.close()

    # t_eddy = tlp[n1:n2];
    t_eddy = [datetime.fromtimestamp(t) for t in tlp[n1:n2]]

    # save relevant arrays
    if not os.path.exists(savedir + estr):
        os.mkdir(savedir + estr)
    pickle.dump(t_eddy, open((savedir + estr + '/t_eddy.p'), 'wb'))
    pickle.dump(bp_anom_eddy, open((savedir + estr + '/bp_eddy.p'), 'wb'))
    pickle.dump(bp_bc_eddy, open((savedir + estr + '/bc_eddy.p'), 'wb'))
    pickle.dump(bp_ssh_eddy, open((savedir + estr + '/ssh_eddy.p'), 'wb'))
    pickle.dump(ctrs_ll, open((savedir + estr + '/ctrs_ll.p'), 'wb'))
    pickle.dump(ctrs_xy, open((savedir + estr + '/ctrs_xy.p'), 'wb'))

# PLOTTING
# eddy track map plot
fig0 = plt.figure(figsize=(11,8.5))

ax1 = fig0.add_subplot(131)
bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
ax1.plot(ctrs_ll[0,0],ctrs_ll[0,1],'ob',label='start')
ax1.plot(ctrs_ll[:,0],ctrs_ll[:,1],'r',label='path')
ax1.plot(ctrs_ll[-1,0],ctrs_ll[-1,1],'xb',label='end')
ax1.legend()
ax1.set_title(estr + ' ' + datetime.fromtimestamp(tlp[n1]).strftime('%m/%d/%y') +
    '-' + datetime.fromtimestamp(tlp[n2]).strftime('%m/%d/%y'))
ax1.grid(True)
ax1.tick_params(labelsize=fs-2)
plt.xticks(rotation=45)

if not os.path.exists(savedir + estr):
    os.mkdir(savedir + estr)

# plt.show()
plt.tight_layout()
plt.savefig(savedir + estr + '/' + estr + '_track')
plt.savefig(savedir + estr + '/' + estr + '_track.eps')
plt.close()

# panel plot of time series at points along track through entire t
for ii in range(25):
    fig2 = plt.figure(figsize=(11,8.5))
    ax2 = fig2.add_subplot(4,3,5)
    iin = int(ii*(n2-n1)/25)
    iilo = ctrs_xy[iin,0]
    iila = ctrs_xy[iin,1]
    ax2.plot(t_eddy,bp_bc_anom[n1:n2,iila,iilo]/100,lw=lw,label='steric')
    ax2.plot(t_eddy,bp_ssh_anom[n1:n2,iila,iilo]/100,lw=lw,label='eustatic')
    ax2.plot(t_eddy,bp_anom2[n1:n2,iila,iilo]/100,lw=lw,color='r',label='total')

    ax2.legend()
    ax2.axhline(c='k',lw=1)
    ax2.axvline(t_eddy[iin],c='k',lw=1)
    myFmt = mdates.DateFormatter('%m/%d')
    ax2.xaxis.set_major_formatter(myFmt)
    plt.xticks(rotation=45, ha='right')
    ax2.set_ylim(-15,15)
    ax2.set_title(str(int(bath[iila,iilo])) + ' m', fontsize=fs)
    ax2.set_ylabel('P (cm)', fontsize=fs)
    ax2.grid(True)
    ax2.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + estr + '/panels'):
        os.mkdir(savedir + estr + '/panels')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + estr + '/panels/' + estr + '_panel_' + str(ii).zfill(2))
    plt.savefig(savedir + estr + '/panels/' + estr + '_panel_' + str(ii).zfill(2) +'.eps')
    plt.close()

# panel plot of differences along track through entire t
for ii in range(25):
    fig2 = plt.figure(figsize=(11,8.5))
    ax2 = fig2.add_subplot(4,3,6)
    iin = int(ii*(n2-n1)/25)
    iiref = int(9*(n2-n1)/25)
    iilo = ctrs_xy[iin,0]
    iila = ctrs_xy[iin,1]
    ax2.plot(t_eddy,(bp_bc_anom[n1:n2,iila,iilo]-bp_bc_anom[n1:n2,ctrs_xy[iiref,1],ctrs_xy[iiref,0]])/100,lw=lw,label='steric')
    ax2.plot(t_eddy,(bp_ssh_anom[n1:n2,iila,iilo]-bp_ssh_anom[n1:n2,ctrs_xy[iiref,1],ctrs_xy[iiref,0]])/100,lw=lw,label='eustatic')
    ax2.plot(t_eddy,(bp_anom2[n1:n2,iila,iilo]-bp_anom2[n1:n2,ctrs_xy[iiref,1],ctrs_xy[iiref,0]])/100,color='r',lw=lw,label='total')

    ax2.legend()
    ax2.axhline(c='k',lw=1)
    ax2.axvline(t_eddy[iin],c='k',lw=1)
    myFmt = mdates.DateFormatter('%m/%d')
    ax2.xaxis.set_major_formatter(myFmt)
    plt.xticks(rotation=45, ha='right')
    ax2.set_ylim(-15,15)
    ax2.set_title(str(int(bath[iila,iilo])) + ' m', fontsize=fs)
    ax2.set_ylabel('\DeltaP (cm)', fontsize=fs)
    ax2.grid(True)
    ax2.tick_params(labelsize=fs-2)

    if not os.path.exists(savedir + estr + '/diff_panels'):
        os.mkdir(savedir + estr + '/diff_panels')

    # plt.show()
    plt.tight_layout()
    plt.savefig(savedir + estr + '/diff_panels/' + estr + '_diffpanel' + str(ii).zfill(2))
    plt.savefig(savedir + estr + '/diff_panels/' + estr + '_diffpanel' + str(ii).zfill(2) + '.eps')
    plt.close()
