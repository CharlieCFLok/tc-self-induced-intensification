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

			ncin_xtime = ncdatain.variables['ocean_time']						# netcdf delta time
			ncin_xlon = ncdatain.variables['lon_rho'][:]						# longitude on mass point [deg]
			ncin_xlat = ncdatain.variables['lat_rho'][:]						# latitude on mass point [deg]

			ncin_temp = ncdatain.variables['temp'][:]						# potential temperature on mass point [C]

			######

			tmp_tstripe = int(1800/(ncin_xtime[1]-ncin_xtime[0]))

			if 'xtime' in tmp_data:
				tmp_data['xtime'] = np.concatenate((tmp_data['xtime'], ncin_xtime[tmp_tstripe::tmp_tstripe]), axis=None)
				tmp_data['temp'] = np.concatenate((tmp_data['temp'], ncin_temp[tmp_tstripe::tmp_tstripe]), axis=0)
			else:
				tmp_data['xtime'] = ncin_xtime[::tmp_tstripe]
				tmp_data['temp'] = ncin_temp[::tmp_tstripe]

		######

		ncin_ntime, ncin_ndepth, ncin_nxlat, ncin_nxlon = tmp_data['temp'].shape

		######

		npzdata_depth = np.load('data/%s_%s/output_proc/ocean_depth.npz' % (tmp_tc, tmp_src, ))
		ncin_depth = npzdata_depth['depth'][:]

		######

		ncin_mld = np.zeros([ncin_ntime, ncin_nxlat, ncin_nxlon])
		ncin_ismld = np.zeros_like(ncin_depth, dtype='bool')

		for tmp_tidx in range(ncin_ntime):
			tmp_mld = tmp_data['temp'][tmp_tidx]-tmp_data['temp'][tmp_tidx, -1]
			tmp_mld[tmp_mld < -0.8] = 9999
			tmp_mldidx = np.argmin(tmp_mld, axis=0)

			###

			for tmp_xlatidx in range(ncin_nxlat):
				for tmp_xlonidx in range(ncin_nxlon):
					ncin_mld[tmp_tidx, tmp_xlatidx, tmp_xlonidx] = ncin_depth[tmp_tidx, tmp_mldidx[tmp_xlatidx, tmp_xlonidx], tmp_xlatidx, tmp_xlonidx]
					ncin_mld[tmp_tidx, tmp_xlatidx, tmp_xlonidx] = min(ncin_mld[tmp_tidx, tmp_xlatidx, tmp_xlonidx], -10)
					ncin_ismld[tmp_tidx, :, tmp_xlatidx, tmp_xlonidx] = (ncin_depth[tmp_tidx, :, tmp_xlatidx, tmp_xlonidx] >= ncin_mld[tmp_tidx, tmp_xlatidx, tmp_xlonidx])

		######

		np.savez('data/%s_%s/output_proc/ocean_mld.npz' % (tmp_tc, tmp_src, ), mld=ncin_mld, ismld=ncin_ismld, )
