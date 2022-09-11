"""
bp_extractor.py
Generate bottom pressure structures using history files on perigee.
Goal is to do so efficiently by only pulling the bare minimum data from
the netCDFs into working memory
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import xarray as xr
import cmocean
import gsw
import pandas as pd
import matplotlib.pyplot as plt
import pickle

sys.path.append(os.path.abspath('util'))
import Lfun
import zrfun
import zfun

# config
testbatch = False # True will only process a few history files, stored locally
if testbatch:
    etag = '_test'
else:
    etag = ''
cutout = True # True will subsample from user-defined lat/lon limits

out_dir = '../LO_output/'
ostr = 'allinone/' # arbitrary label for separating runs
in_dir = '/data1/parker/LO_roms/cas6_v0_live/' # location on perigee
# in_dir = '../LO_data/cas6_v0_live/' # location locally
dstr = 'f' # naming convention for directories
fstr = 'ocean_his_' # naming convention for history files
ti = datetime.strptime('2020.12.15', '%Y.%m.%d')
tf = datetime.strptime('2021.12.14', '%Y.%m.%d')
# tf = datetime.strptime('2022.06.30', '%Y.%m.%d')
tag = 'LiveOcean'
n_layer = 0 # bottom layer

# gridname = 'cascadia1'
# tag = 'base'
# ex_name = 'lobio5'
# list_type = 'low_passed'
#
# Ldir = Lfun.Lstart(gridname, tag)
# Ldir['ex_name'] = ex_name
# Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']
# tag = Ldir['gtagex']
#
# if Ldir['lo_env'] == 'pm_mac':
#     dt0 = datetime(2017,8,5)
#     dt1 = datetime(2017,8,9)
# elif Ldir['lo_env'] == 'pm_fjord':
#     dt0 = datetime(2017,1,1)
#     dt1 = datetime(2017,12,31)
# #
# if list_type == 'low_passed':
#     fn_list = []
#     dt = dt0
#     while dt <= dt1:
#         date_string = dt.strftime(format='%Y.%m.%d')
#         Ldir['date_string'] = date_string
#         f_string = 'f' + Ldir['date_string']
#         in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
#         if 'low_passed.nc' in os.listdir(in_dir):
#             fn_list.append(in_dir + 'low_passed.nc')
#         else:
#             print('Missing file for ' + date_string)
#         dt = dt + timedelta(days=1)
#

#### end setup ####

# get grid info first (things that won't change between files)
file0 = in_dir + dstr + datetime.strftime(ti, '%Y.%m.%d') + '/' + fstr + '0001.nc'
ds0 = nc.Dataset(file0) # opens netCDF file
if False: # set to True if you want to print netCDF info
    print('\nGRID INFO')
    zfun.ncd(ds0)
lon = ds0['lon_rho'][:,:]
lat = ds0['lat_rho'][:,:]
bath = ds0['h'][:,:]
if cutout: # take cutout from full model domain, if desired
    # full range is: -129.9798<lon<-122.018, 42.0067<lat<52.0099
    minlon = -127; maxlon = -124; minlat = 44; maxlat = 48
    ilo1 = np.argmin(np.abs(lon[0,:] - minlon)); ilo2 = np.argmin(np.abs(lon[0,:] - maxlon))
    ila1 = np.argmin(np.abs(lat[:,0] - minlat)); ila2 = np.argmin(np.abs(lat[:,0] - maxlat))
    lon = lon[ila1:ila2,ilo1:ilo2]
    lat = lat[ila1:ila2,ilo1:ilo2]
    bath = bath[ila1:ila2,ilo1:ilo2]
ds0.close()

# make file lists
file_list = []
if testbatch: # only runs on 2 day subset
    delta = (ti + timedelta(days=2)) - ti
else: # runs on entire date range defined in CONFIG
    delta = (tf + timedelta(days=1)) - ti
drange = range(0, delta.days)
for ii in drange:  # building file names
    d_ii = in_dir + dstr + datetime.strftime(ti + timedelta(days=ii), '%Y.%m.%d') + '/'
    for jj in range(1,25):
        f_jj = fstr + str(jj).zfill(4) + '.nc'
        file_list.append(d_ii + f_jj)
nf = len(file_list)

# make some things
fn = file_list[0]
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = ds['h'][:]
z = zrfun.get_z(h, 0*h, S, only_rho=True)
z0 = z[n_layer,:,:].squeeze()
ds.close()

# prepare a directory for results
outdir0 = out_dir + ostr
Lfun.make_dir(outdir0, clean=False)
ncoutdir = outdir0 + 'pickles/'
Lfun.make_dir(ncoutdir, clean=False)
# output file
out_name = 'bp_' + tag + etag + '.nc'
out_fn = ncoutdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# make bottom pressure
bp_tot = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
bp_tot2 = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
bp_bc = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
bp_ssh = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
ssh_tot = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
t_arr = np.zeros(nf)

"""
# get time axis
ot = ds.ocean_time.values # an array with dtype='datetime64[ns]'
dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
dt = dti.to_pydatetime() # an array of datetimes
"""

# set constants
pad = 36 # this trims the ends after the low pass so there are no nan's
g = 9.81 # gravity [m s-2]
for tt in range(nf):
    # if tt==59 or tt==60 or tt==63 or tt==70: # something is wrong with these files (at least on local)
    #     bp_arr[tt,:,:] = bp_arr[tt-1,:,:] # populate with previous day
    #     t_arr[tt] = t_arr[tt-1]
    #     continue
    filex = file_list[tt]
    dsx = nc.Dataset(filex)
    if np.mod(tt,24)==0: # print update to console for every day
        print('tt = ' + str(tt) + '/' + str(nf) + ' ' + filex[len(filex)-28:len(filex)])
        sys.stdout.flush()
    t_arr[tt] = dsx['ocean_time'][0].squeeze()
    if cutout:
        zeta = dsx['zeta'][0, ila1:ila2, ilo1:ilo2].squeeze()
        rho = dsx['rho'][0, :, ila1:ila2, ilo1:ilo2].squeeze() + 1000.
        salt = dsx['salt'][0, :, ila1:ila2, ilo1:ilo2].squeeze()
        temp = dsx['temp'][0, :, ila1:ila2, ilo1:ilo2].squeeze()
        z_rho,z_w = zrfun.get_z(G['h'][ila1:ila2, ilo1:ilo2], zeta, S)
    else:
        zeta = dsx['zeta'][0,:,:].squeeze()
        rho = dsx['rho'][0,:,:,:].squeeze() + 1000.
        salt = dsx['salt'][0,:,:,:].squeeze()
        temp = dsx['temp'][0,:,:,:].squeeze()
        z_rho,z_w = zrfun.get_z(G['h'], zeta, S)
    if tt == 0:
        bp_arr = (0*zeta) * np.ones((nf,1,1))
        DA = G['DX'] * G['DY']
        # DAm = np.ma.masked_where(zeta.mask, DA)
    NT = 1
    N = 30

    Z = z_rho - z_w[-1] # adjust so free surface is at 0
    DZ = np.diff(z_w, axis=0)
    bp_tot[tt,:,:] = (g * rho * DZ).sum(axis=0)
    ssh_tot[tt,:,:] = zeta

    # Equation of state calculations
    p = gsw.p_from_z(Z, lat)
    SA = gsw.SA_from_SP(salt, p, lon, lat) # absolute salinity
    CT = gsw.CT_from_pt(SA, temp) # conservative temperature
    if True:
        rho = gsw.rho(SA, CT, p)
        # This is denser than ROMS rho by 0.037 [kg m-3] at the bottom and 0.0046 [kg m-3]
        # (annual averages), and it is the full density, not density minus 1000.
        # There was no visual difference between the pressure time series.

    # calculate the baroclinic pressure
    p = np.flip(np.cumsum(np.flip(g * rho * DZ, axis=0), axis=0), axis=0)

    # calculate the pressure due to SSH
    p0 = g * 1025 * zeta

    # # separate contributions of salt and temp to density
    # rho_only_salt = gsw.rho(SA, CT - CT, p) - Rho.reshape((1,N))
    # rho_only_temp = gsw.rho(SA - SA, CT, p) - Rho.reshape((1,N))

    # make the full pressure anomaly
    bp_tot2[tt,:,:] = p[0,:,:] + p0
    bp_bc[tt,:,:] = p[0,:,:]
    bp_ssh[tt,:,:] = p0

    dsx.close()

"""
# convert to anomalies
sshm = np.mean(ssh_tot, axis=0)
ssh_anom = ssh_tot - sshm
bpm = np.mean(bp_tot, axis=0)
bp_anom = bp_tot - bpm
bpm = np.mean(bp_tot2, axis=0)
bp_anom2 = bp_tot2 - bpm
bpm = np.mean(bp_bc, axis=0)
bp_bc_anom = bp_bc - bpm
bpm = np.mean(bp_ssh, axis=0)
bp_ssh_anom = bp_ssh - bpm

# annual means
Rho = np.mean(rho, axis=0)
Salt = np.mean(SA, axis=0)
Temp = np.mean(CT, axis=0)
P = np.mean(p, axis=0)

# anomalies from the annual mean
salt_a = SA - Salt.reshape((1,N))
temp_a = CT - Temp.reshape((1,N))
rho_a = rho - Rho.reshape((1,N))
p_a = p - P.reshape((1,N))

# low pass filtered version
zetalp = zfun.lowpass(zeta, f='godin')[pad:-pad:24]
zetalp = etalp - np.mean(zetalp) # remove mean SSH
rholp = zfun.lowpass(rho, f='godin')[pad:-pad:24, :]
saltlp = zfun.lowpass(SA, f='godin')[pad:-pad:24, :]
templp = zfun.lowpass(CT, f='godin')[pad:-pad:24, :]

# also make associated time vectors
tlp = dt[pad:-pad:24]
tlpf = dt[24:-24:24] # "full" version used for pcolormesh plots
NTlp = len(etalp)
"""

# SAVING
pickle.dump(t_arr, open((ncoutdir + 't_arr.p'), 'wb'))
pickle.dump(bp_tot, open((ncoutdir + 'bp_tot.p'), 'wb'))
# pickle.dump(bp_anom, open((ncoutdir + 'bp_anom.p'), 'wb'))
pickle.dump(bp_tot2, open((ncoutdir + 'bp_tot2.p'), 'wb'))
# pickle.dump(bp_anom2, open((ncoutdir + 'bp_anom2.p'), 'wb'))
pickle.dump(bp_bc, open((ncoutdir + 'bp_bc.p'), 'wb'))
# pickle.dump(bp_bc_anom, open((ncoutdir + 'bp_bc_anom.p'), 'wb'))
pickle.dump(bp_ssh, open((ncoutdir + 'bp_ssh.p'), 'wb'))
# pickle.dump(bp_ssh_anom, open((ncoutdir + 'bp_ssh_anom.p'), 'wb'))
pickle.dump(lat, open((ncoutdir + 'lat.p'), 'wb'))
pickle.dump(lon, open((ncoutdir + 'lon.p'), 'wb'))
