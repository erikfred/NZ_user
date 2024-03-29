"""
map_fields.py
Generate mapview plots of desired fields.
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import cmocean
import matplotlib.pyplot as plt
import pickle
from scipy import signal

from lo_tools import Lfun, zfun, zrfun

g = 9.81
topdir = '../LO_output/allinone/'
loadir = topdir + 'pickles_2019-20/'
outdir = '../LO_output/mapview/'
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

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.balance # formerly thermal
cmap2 = cmocean.cm.haline

plt.close('all')

# # RMS OF BOTTOM PRESSURE COMPONENT ANOMALIES
# levels=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
# # total anomaly
# fig1 = plt.figure(figsize=(8,8))
# ax1 = fig1.add_subplot(111)
# cs = ax1.contourf(lon, lat, np.std(bp_anom2, axis=0)/100, levels=levels, cmap=cmap2, extend='max')
# bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
# ax1.axis('square')
# # ax1.set_xlim((minlon,maxlon))
# # ax1.set_ylim((minlat, maxlat))
# ax1.grid(True)
# ax1.tick_params(labelsize=fs-2)
# plt.xticks(rotation=45)
#
# cb = fig1.colorbar(cs)
# cb.set_label('RMS (cm)', fontsize=fs-2)
# cb.ax.tick_params(labelsize=fs-2)
#
# # plt.show()
# if not os.path.exists(savedir + 'bp_anom'):
#     os.mkdir(savedir + 'bp_anom')
# plt.savefig(savedir + 'bp_anom/bp_anom_rms')
# plt.savefig(savedir + 'bp_anom/bp_anom_rms.eps')
# plt.close()
#
# # baroclinic component
# fig2 = plt.figure(figsize=(8,8))
# ax2 = fig2.add_subplot(111)
# cs = ax2.contourf(lon, lat, np.std(bp_bc_anom, axis=0)/100, levels=levels, cmap=cmap2, extend='max')
# bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
# ax2.axis('square')
# ax2.grid(True)
# ax2.tick_params(labelsize=fs-2)
# plt.xticks(rotation=45)
#
# cb = fig2.colorbar(cs)
# cb.set_label('RMS (cm)', fontsize=fs-2)
# cb.ax.tick_params(labelsize=fs-2)
#
# # plt.show()
# if not os.path.exists(savedir + 'bp_bc'):
#     os.mkdir(savedir + 'bp_bc')
# plt.savefig(savedir + 'bp_bc/bp_bc_rms')
# plt.savefig(savedir + 'bp_bc/bp_bc_rms.eps')
# plt.close()
#
# # ssh component
# fig3 = plt.figure(figsize=(8,8))
# ax3 = fig3.add_subplot(111)
# cs = ax3.contourf(lon, lat, np.std(bp_ssh_anom, axis=0)/100, levels=levels, cmap=cmap2, extend='max')
# bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
# ax3.axis('square')
# ax3.grid(True)
# ax3.tick_params(labelsize=fs-2)
# plt.xticks(rotation=45)
#
# cb = fig3.colorbar(cs)
# cb.set_label('RMS (cm)', fontsize=fs-2)
# cb.ax.tick_params(labelsize=fs-2)
#
# # plt.show()
# if not os.path.exists(savedir + 'bp_ssh'):
#     os.mkdir(savedir + 'bp_ssh')
# plt.savefig(savedir + 'bp_ssh/bp_ssh_rms')
# plt.savefig(savedir + 'bp_ssh/bp_ssh_rms.eps')
# plt.close()

# # DAILY BOTTOM PRESSURE COMPONENT ANOMALIES
# for nn in range(len(bp_anom2[:,0,0])): # for each day
#     di = datetime.fromtimestamp(tlp[nn])
#     # total anomaly
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.pcolormesh(lon, lat, bp_anom2[nn,:,:]/100, cmap=cmap, vmin=-15, vmax=15)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     # ax1.set_xlim((minlon,maxlon))
#     # ax1.set_ylim((minlat, maxlat))
#     ax1.grid(True)
#     ax1.set_title('BPA ' + di.strftime('%m/%d/%y') + ' (n = ' + str(nn) + ')')
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('Pressure (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_anom'):
#         os.mkdir(savedir + 'bp_anom')
#     plt.savefig(savedir + 'bp_anom/bp_anom_' + str(nn).zfill(4))
#     plt.close()
#
#     # baroclinic component
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.pcolormesh(lon, lat, bp_bc_anom[nn,:,:]/100, cmap=cmap, vmin=-15, vmax=15)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#     ax2.set_title('P_BC ' + di.strftime('%m/%d/%y') + ' (n = ' + str(nn) + ')')
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('Pressure (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_bc'):
#         os.mkdir(savedir + 'bp_bc')
#     plt.savefig(savedir + 'bp_bc/bp_bc_' + str(nn).zfill(4))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.pcolormesh(lon, lat, bp_ssh_anom[nn,:,:]/100, cmap=cmap, vmin=-15, vmax=15)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#     ax3.set_title('P_SSH ' + di.strftime('%m/%d/%y') + ' (n = ' + str(nn) + ')')
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('Pressure (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_ssh'):
#         os.mkdir(savedir + 'bp_ssh')
#     plt.savefig(savedir + 'bp_ssh/bp_ssh_' + str(nn).zfill(4))
#     plt.close()

# # CORRELATION COEFFICIENTS RELATIVE TO REFERENCE POINTS
# # strike-parallel line
# c_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_bc = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_ssh = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# isb = np.argwhere((bath<302) & (bath>298))
# for nn in range(int(np.around(len(isb)/5))):
#     for mm in range(len(bp_anom2[0,0,:])): # for every grid point
#         for kk in range(len(bp_anom2[0,:,0])): # for every grid point
#             # full anomaly
#             p1 = bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_anom2[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_anom[kk,mm] = np.correlate(p1, p2)
#             # baroclinic component
#             p1 = bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_bc_anom[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_bc[kk,mm] = np.correlate(p1, p2)
#             # ssh component
#             p1 = bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_ssh_anom[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_ssh[kk,mm] = np.correlate(p1, p2)
#
#     # total anomaly
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.pcolormesh(lon, lat, c_anom, cmap=cmap, vmin=-1, vmax=1)
#     ax1.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     ax1.grid(True)
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_anom'):
#         os.mkdir(savedir + 'bp_anom')
#     plt.savefig(savedir + 'bp_anom/xcorr_strikepar_' + str(nn))
#     plt.close()
#
#     # baroclinic component
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.pcolormesh(lon, lat, c_bc, cmap=cmap, vmin=-1, vmax=1)
#     ax2.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_bc'):
#         os.mkdir(savedir + 'bp_bc')
#     plt.savefig(savedir + 'bp_bc/xcorr_strikepar_' + str(nn))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.pcolormesh(lon, lat, c_ssh, cmap=cmap, vmin=-1, vmax=1)
#     ax3.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_ssh'):
#         os.mkdir(savedir + 'bp_ssh')
#     plt.savefig(savedir + 'bp_ssh/xcorr_strikepar_' + str(nn))
#     plt.close()
#
# # strike-perpendicular line
# for nn in range(20):
#     for mm in range(len(bp_anom2[0,0,:])): # for every grid point
#         for kk in range(len(bp_anom2[0,:,0])): # for every grid point
#             # full anomaly
#             p1 = bp_anom2[:,30,nn*10]
#             p2 = bp_anom2[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_anom[kk,mm] = np.correlate(p1, p2)
#             # baroclinic component
#             p1 = bp_bc_anom[:,30,nn*10]
#             p2 = bp_bc_anom[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_bc[kk,mm] = np.correlate(p1, p2)
#             # ssh component
#             p1 = bp_ssh_anom[:,30,nn*10]
#             p2 = bp_ssh_anom[:,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_ssh[kk,mm] = np.correlate(p1, p2)
#
#     # total anomaly
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.pcolormesh(lon, lat, c_anom, cmap=cmap, vmin=-1, vmax=1)
#     ax1.plot(lon[30,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     ax1.grid(True)
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_anom'):
#         os.mkdir(savedir + 'bp_anom')
#     plt.savefig(savedir + 'bp_anom/xcorr_strikeperp_' + str(nn))
#     plt.close()
#
#     # baroclinic component
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.pcolormesh(lon, lat, c_bc, cmap=cmap, vmin=-1, vmax=1)
#     ax2.plot(lon[30,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_bc'):
#         os.mkdir(savedir + 'bp_bc')
#     plt.savefig(savedir + 'bp_bc/xcorr_strikeperp_' + str(nn))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.pcolormesh(lon, lat, c_ssh, cmap=cmap, vmin=-1, vmax=1)
#     ax3.plot(lon[30,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_ssh'):
#         os.mkdir(savedir + 'bp_ssh')
#     plt.savefig(savedir + 'bp_ssh/xcorr_strikeperp_' + str(nn))
#     plt.close()

# # SAME LINES, NOW RMS AS METRIC
# # strike-parallel line
# c_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_bc = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_ssh = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# isb = np.argwhere((bath<502) & (bath>498))
# for nn in range(int(np.around(len(isb)/5))):
#     for mm in range(len(bp_anom2[0,0,:])): # for every grid point
#         for kk in range(len(bp_anom2[0,:,0])): # for every grid point
#             # full anomaly
#             p1 = bp_anom2[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_anom2[:,kk,mm]
#             c_anom[kk,mm] = np.std(p1-p2)
#             # baroclinic component
#             p1 = bp_bc_anom[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_bc_anom[:,kk,mm]
#             c_bc[kk,mm] = np.std(p1-p2)
#             # ssh component
#             p1 = bp_ssh_anom[:,isb[nn*5][0],isb[nn*5][1]]
#             p2 = bp_ssh_anom[:,kk,mm]
#             c_ssh[kk,mm] = np.std(p1-p2)
#
#     # total anomaly
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.contourf(lon, lat, c_anom/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax1.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     ax1.grid(True)
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_anom'):
#         os.mkdir(savedir + 'bp_anom')
#     plt.savefig(savedir + 'bp_anom/rms_strikepar500_' + str(nn))
#     plt.close()
#
#     # baroclinic component
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.contourf(lon, lat, c_bc/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax2.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_bc'):
#         os.mkdir(savedir + 'bp_bc')
#     plt.savefig(savedir + 'bp_bc/rms_strikepar500_' + str(nn))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.contourf(lon, lat, c_ssh/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax3.plot(lon[isb[nn*5][0],isb[nn*5][1]], lat[isb[nn*5][0],isb[nn*5][1]], '*k', markersize=mk)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_ssh'):
#         os.mkdir(savedir + 'bp_ssh')
#     plt.savefig(savedir + 'bp_ssh/rms_strikepar500_' + str(nn))
#     plt.close()
#
# # strike-perpendicular line
# for nn in range(5,24):
#     for mm in range(len(bp_anom2[0,0,:])): # for every grid point
#         for kk in range(len(bp_anom2[0,:,0])): # for every grid point
#             # full anomaly
#             p1 = bp_anom2[:,100,nn*10]
#             p2 = bp_anom2[:,kk,mm]
#             c_anom[kk,mm] = np.std(p1-p2)
#             # baroclinic component
#             p1 = bp_bc_anom[:,100,nn*10]
#             p2 = bp_bc_anom[:,kk,mm]
#             c_bc[kk,mm] = np.std(p1-p2)
#             # ssh component
#             p1 = bp_ssh_anom[:,100,nn*10]
#             p2 = bp_ssh_anom[:,kk,mm]
#             c_ssh[kk,mm] = np.std(p1-p2)
#
#     # total anomaly
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.contourf(lon, lat, c_anom/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax1.plot(lon[100,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     ax1.grid(True)
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_anom'):
#         os.mkdir(savedir + 'bp_anom')
#     plt.savefig(savedir + 'bp_anom/rms_strikeperp_' + str(nn))
#     plt.close()
#
#     # baroclinic component
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.contourf(lon, lat, c_bc/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax2.plot(lon[100,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_bc'):
#         os.mkdir(savedir + 'bp_bc')
#     plt.savefig(savedir + 'bp_bc/rms_strikeperp_' + str(nn))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.contourf(lon, lat, c_ssh/100, levels=[0, 0.5, 1, 1.5, 2, 2.5], cmap=cmap2)
#     ax3.plot(lon[100,nn*10], lat[30,nn*10], '*k', markersize=mk)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('RMS (cm)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir + 'bp_ssh'):
#         os.mkdir(savedir + 'bp_ssh')
#     plt.savefig(savedir + 'bp_ssh/rms_strikeperp_' + str(nn))
#     plt.close()

# CORRELATION COEFFICIENTS BETWEEN COMPONENTS
c_ssh_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
c_bc_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
c_ssh_bc = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
for mm in range(len(bp_anom2[0,0,:])): # for every grid point
    for kk in range(len(bp_anom2[0,:,0])): # for every grid point
        # ssh / anom
        p1 = bp_ssh_anom[:,kk,mm]
        p2 = bp_anom2[:,kk,mm]
        p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
        p2 = (p2 - np.mean(p2)) / (np.std(p2))
        c_ssh_anom[kk,mm] = np.correlate(p1, p2)
        # baroclinic / anom
        p1 = bp_bc_anom[:,kk,mm]
        p2 = bp_anom2[:,kk,mm]
        p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
        p2 = (p2 - np.mean(p2)) / (np.std(p2))
        c_bc_anom[kk,mm] = np.correlate(p1, p2)
        # ssh / baroclinic
        p1 = bp_ssh_anom[:,kk,mm]
        p2 = bp_bc_anom[:,kk,mm]
        p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
        p2 = (p2 - np.mean(p2)) / (np.std(p2))
        c_ssh_bc[kk,mm] = np.correlate(p1, p2)

levels=[-1.0, -0.9, -0.8, -0.7, -0.6, 0.6, 0.7, 0.8, 0.9, 1.0]
# ssh / anom
mask = np.absolute(c_bc_anom) > 0.2
C = np.ma.masked_array(c_ssh_anom, mask=mask)
fig1 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(111)
# cs = ax1.pcolormesh(lon, lat, c_ssh_anom, cmap=cmap, vmin=-1, vmax=1)
cs = ax1.contourf(lon, lat, c_ssh_anom, levels=levels, cmap=cmap)
bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')

ax1.axis('square')
ax1.grid(True)
ax1.tick_params(labelsize=fs-2)
plt.xticks(rotation=45)

cb = fig1.colorbar(cs)
cb.set_label('Cross Correlation', fontsize=fs-2)
cb.ax.tick_params(labelsize=fs-2)

# plt.show()
if not os.path.exists(savedir + 'cross_field'):
    os.mkdir(savedir + 'cross_field')
plt.savefig(savedir + 'cross_field/xcorr_ssh_anom')
plt.savefig(savedir + 'cross_field/xcorr_ssh_anom.eps')
# cs.remove()
# cs = ax1.pcolormesh(lon, lat, C, cmap=cmap, vmin=-1, vmax=1)
# ax1.grid(True)
# plt.savefig(savedir + 'cross_field/xcorr_ssh_anom_masked')
plt.close()

# bc / anom
mask = np.absolute(c_ssh_anom) > 0.2
C = np.ma.masked_array(c_bc_anom, mask=mask)
fig2 = plt.figure(figsize=(8,8))
ax2 = fig2.add_subplot(111)
# cs = ax2.pcolormesh(lon, lat, c_bc_anom, cmap=cmap, vmin=-1, vmax=1)
cs = ax2.contourf(lon, lat, c_bc_anom, levels=levels, cmap=cmap)
bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')

ax2.axis('square')
ax2.grid(True)
ax2.tick_params(labelsize=fs-2)
plt.xticks(rotation=45)

cb = fig2.colorbar(cs)
cb.set_label('Cross Correlation', fontsize=fs-2)
cb.ax.tick_params(labelsize=fs-2)

# plt.show()
if not os.path.exists(savedir + 'cross_field'):
    os.mkdir(savedir + 'cross_field')
plt.savefig(savedir + 'cross_field/xcorr_bc_anom')
plt.savefig(savedir + 'cross_field/xcorr_bc_anom.eps')
# cs.remove()
# cs = ax2.pcolormesh(lon, lat, C, cmap=cmap, vmin=-1, vmax=1)
# ax2.grid(True)
# plt.savefig(savedir + 'cross_field/xcorr_bc_anom_masked')
plt.close()

# ssh / bc
fig3 = plt.figure(figsize=(8,8))
ax3 = fig3.add_subplot(111)
# cs = ax3.pcolormesh(lon, lat, c_ssh_bc, cmap=cmap, vmin=-1, vmax=1)
cs = ax3.contourf(lon, lat, c_ssh_bc, levels=levels, cmap=cmap)
bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')

ax3.axis('square')
ax3.grid(True)
ax3.tick_params(labelsize=fs-2)
plt.xticks(rotation=45)

cb = fig3.colorbar(cs)
cb.set_label('Cross Correlation', fontsize=fs)
cb.ax.tick_params(labelsize=fs)

# plt.show()
if not os.path.exists(savedir + 'cross_field'):
    os.mkdir(savedir + 'cross_field')
plt.savefig(savedir + 'cross_field/xcorr_ssh_bc')
plt.savefig(savedir + 'cross_field/xcorr_ssh_bc.eps')
plt.close()

# # AS ABOVE, BUT WINDOWED RATHER THAN FULL-YEAR
# win = 60 # play around with this number
# savedir2 = savedir + 'cross_field/'
# c_ssh_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_bc_anom = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# c_ssh_bc = np.empty((np.shape(lon)[0], np.shape(lon)[1]),dtype=np.float32);
# for nn in range(len(bp_anom2[:,0,0])-win):
#     di = datetime.fromtimestamp(tlp[nn])
#     df = datetime.fromtimestamp(tlp[nn+win])
#     for mm in range(len(bp_anom2[0,0,:])): # for every grid point
#         for kk in range(len(bp_anom2[0,:,0])): # for every grid point
#             # ssh / anom
#             p1 = bp_ssh_anom[nn:nn+win,kk,mm]
#             p2 = bp_anom2[nn:nn+win,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_ssh_anom[kk,mm] = np.correlate(p1, p2)
#             # baroclinic / anom
#             p1 = bp_bc_anom[nn:nn+win,kk,mm]
#             p2 = bp_anom2[nn:nn+win,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_bc_anom[kk,mm] = np.correlate(p1, p2)
#             # ssh / baroclinic
#             p1 = bp_ssh_anom[nn:nn+win,kk,mm]
#             p2 = bp_bc_anom[nn:nn+win,kk,mm]
#             p1 = (p1 - np.mean(p1)) / (np.std(p1) * len(p1))
#             p2 = (p2 - np.mean(p2)) / (np.std(p2))
#             c_ssh_bc[kk,mm] = np.correlate(p1, p2)
#
#     # ssh / anom
#     mask = np.absolute(c_bc_anom) > 0.2
#     C = np.ma.masked_array(c_ssh_anom, mask=mask)
#     fig1 = plt.figure(figsize=(8,8))
#     ax1 = fig1.add_subplot(111)
#     cs = ax1.pcolormesh(lon, lat, c_ssh_anom, cmap=cmap, vmin=-1, vmax=1)
#     bth = ax1.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax1.axis('square')
#     ax1.grid(True)
#     ax1.set_title('SSH v. ANOM ' + di.strftime('%m/%d') + ' - ' + df.strftime('%m/%d/%Y'))
#
#     cb = fig1.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir2 + str(win) + 'day'):
#         os.mkdir(savedir2 + str(win) + 'day')
#     plt.savefig(savedir2 + str(win) + 'day/xcorr_ssh_anom_' + str(nn).zfill(4))
#     cs.remove()
#     cs = ax1.pcolormesh(lon, lat, C, cmap=cmap, vmin=-1, vmax=1)
#     ax1.grid(True)
#     plt.savefig(savedir2 + str(win) + 'day/xcorr_ssh_anom_masked_' + str(nn).zfill(4))
#     plt.close()
#
#     # baroclinic component
#     mask = np.absolute(c_ssh_anom) > 0.2
#     C = np.ma.masked_array(c_bc_anom, mask=mask)
#     fig2 = plt.figure(figsize=(8,8))
#     ax2 = fig2.add_subplot(111)
#     cs = ax2.pcolormesh(lon, lat, c_bc_anom, cmap=cmap, vmin=-1, vmax=1)
#     bth = ax2.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax2.axis('square')
#     ax2.grid(True)
#     ax2.set_title('BC v. ANOM ' + di.strftime('%m/%d') + ' - ' + df.strftime('%m/%d/%Y'))
#
#     cb = fig2.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir2 + str(win) + 'day'):
#         os.mkdir(savedir2 + str(win) + 'day')
#     plt.savefig(savedir2 + str(win) + 'day/xcorr_bc_anom_' + str(nn).zfill(4))
#     cs.remove()
#     cs = ax2.pcolormesh(lon, lat, C, cmap=cmap, vmin=-1, vmax=1)
#     ax2.grid(True)
#     plt.savefig(savedir2 + str(win) + 'day/xcorr_bc_anom_masked_' + str(nn).zfill(4))
#     plt.close()
#
#     # ssh component
#     fig3 = plt.figure(figsize=(8,8))
#     ax3 = fig3.add_subplot(111)
#     cs = ax3.pcolormesh(lon, lat, c_ssh_bc, cmap=cmap, vmin=-1, vmax=1)
#     bth = ax3.contour(lon, lat, bath, [4, 300, 2000], colors='black')
#
#     ax3.axis('square')
#     ax3.grid(True)
#     ax3.set_title('SSH v. BC ' + di.strftime('%m/%d') + ' - ' + df.strftime('%m/%d/%Y'))
#
#     cb = fig3.colorbar(cs)
#     cb.set_label('Cross Correlation', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(savedir2 + str(win) + 'day'):
#         os.mkdir(savedir2 + str(win) + 'day')
#     plt.savefig(savedir2 + str(win) + 'day/xcorr_ssh_bc_' + str(nn).zfill(4))
#     plt.close()
