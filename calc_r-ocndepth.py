#!/home/clok/anaconda3/bin/python3

### system
import sys, os
### netcdf
from netCDF4 import Dataset
#import wrf
### numerical
import numpy as np
#from scipy import interpolate as sp_intpl
### calendar
#import datetime as dt
#import pandas as pd

### plotting
#import matplotlib.pyplot as plt
#import matplotlib.dates as plt_dates
#import cartopy.crs as plt_mcrs
#import cartopy.feature as plt_mfeature

######

plot_tc = ['HATO', 'MANGKHUT', 'NIDA', ]
plot_run = ['ctrl', 'r0rad', 'r0rad-to18', ]

######

for tmp_tc in plot_tc:
	for tmp_src in plot_run:
		tmp_data = {}
		
		for tmp_fnnc in os.listdir('data/%s_%s/output_nc/' % (tmp_tc, tmp_src, )):
			if 'ocean_his_d01' not in tmp_fnnc:
				continue

			print('Opening NetCDF file: data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))
			ncdatain = Dataset('data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))

			###

			ncin_xtime = ncdatain.variables['ocean_time']						# netcdf delta time
			ncin_xlon = ncdatain.variables['lon_rho'][:]						# longitude on mass point [deg]
			ncin_xlat = ncdatain.variables['lat_rho'][:]						# latitude on mass point [deg]

			ncin_zs = ncdatain.variables['s_rho'][:]
			ncin_zs_w = ncdatain.variables['s_w'][:]

			ncin_cs = ncdatain.variables['Cs_r'][:]
			ncin_cs_w = ncdatain.variables['Cs_w'][:]

			ncin_hc = ncdatain.variables['hc'][:]
			ncin_bathy = ncdatain.variables['h'][:]
			ncin_zeta = ncdatain.variables['zeta'][:]

			######

			tmp_tstripe = int(1800/(ncin_xtime[1]-ncin_xtime[0]))

			if 'xtime' in tmp_data:
				tmp_data['xtime'] = np.concatenate((tmp_data['xtime'], ncin_xtime[tmp_tstripe::tmp_tstripe]), axis=None)
				tmp_data['zeta'] = np.concatenate((tmp_data['zeta'], ncin_zeta[tmp_tstripe::tmp_tstripe]), axis=0)
			else:
				tmp_data['xtime'] = ncin_xtime[::tmp_tstripe]
				tmp_data['zeta'] = ncin_zeta[::tmp_tstripe]

		######

		ncin_nxlon = ncin_xlon.shape[1]
		ncin_nxlat = ncin_xlat.shape[0]
		ncin_ndepth = ncin_zs.shape[0]
		ncin_ndepth_w = ncin_zs_w.shape[0]
		ncin_ntime = tmp_data['xtime'].shape[0]

		######

		ncin_depth = np.zeros([ncin_ntime, ncin_ndepth, ncin_nxlat, ncin_nxlon])
		ncin_depth_w = np.zeros([ncin_ntime, ncin_ndepth_w, ncin_nxlat, ncin_nxlon])

		for tmp_tidx in range(ncin_ntime):
			for tmp_zidx in range(ncin_ndepth):
				ncin_depth[tmp_tidx, tmp_zidx, :, :] = tmp_data['zeta'][tmp_tidx, :, :]+(tmp_data['zeta'][tmp_tidx, :, :]+ncin_bathy)*((ncin_hc*ncin_zs[tmp_zidx]+ncin_cs[tmp_zidx]*ncin_bathy)/(ncin_hc+ncin_bathy))
			for tmp_zidx in range(ncin_ndepth_w):
				ncin_depth_w[tmp_tidx, tmp_zidx, :, :] = tmp_data['zeta'][tmp_tidx, :, :]+(tmp_data['zeta'][tmp_tidx, :, :]+ncin_bathy)*((ncin_hc*ncin_zs_w[tmp_zidx]+ncin_cs_w[tmp_zidx]*ncin_bathy)/(ncin_hc+ncin_bathy))

		######

		np.savez('data/%s_%s/output_proc/ocean_depth.npz' % (tmp_tc, tmp_src, ), depth=ncin_depth, depth_w=ncin_depth_w)
