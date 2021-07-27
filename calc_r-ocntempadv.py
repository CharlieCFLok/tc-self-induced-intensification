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

const = {	'a':	6371009,						# [m]
		'g0':	9.81,							# [m s**(-1)]
		}

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

			ncin_ucurr = ncdatain.variables['u_eastward'][:]					# earth rotated u-current on mass point [m s**(-1)]
			ncin_vcurr = ncdatain.variables['v_northward'][:]					# earth rotated u-current on mass point [m s**(-1)]
			ncin_w = ncdatain.variables['w'][:]							# w-current on mass point [m s**(-1)]

			######

			tmp_tstripe = int(1800/(ncin_xtime[1]-ncin_xtime[0]))

			if 'xtime' in tmp_data:
				tmp_data['xtime'] = np.concatenate((tmp_data['xtime'], ncin_xtime[tmp_tstripe::tmp_tstripe]), axis=None)
				tmp_data['temp'] = np.concatenate((tmp_data['temp'], ncin_temp[tmp_tstripe::tmp_tstripe]), axis=0)
				tmp_data['ucurr'] = np.concatenate((tmp_data['ucurr'], ncin_ucurr[tmp_tstripe::tmp_tstripe]), axis=0)
				tmp_data['vcurr'] = np.concatenate((tmp_data['vcurr'], ncin_vcurr[tmp_tstripe::tmp_tstripe]), axis=0)
				tmp_data['w'] = np.concatenate((tmp_data['w'], ncin_w[tmp_tstripe::tmp_tstripe]), axis=0)
			else:
				tmp_data['xtime'] = ncin_xtime[::tmp_tstripe]
				tmp_data['temp'] = ncin_temp[::tmp_tstripe]
				tmp_data['ucurr'] = ncin_ucurr[::tmp_tstripe]
				tmp_data['vcurr'] = ncin_vcurr[::tmp_tstripe]
				tmp_data['w'] = ncin_w[::tmp_tstripe]

		######

		npzdata_depth = np.load('data/%s_%s/output_proc/ocean_depth.npz' % (tmp_tc, tmp_src, ))
		ncin_depth = npzdata_depth['depth'][:]

		######

		tmp_rlon = np.radians(ncin_xlon)
		tmp_rlat = np.radians(ncin_xlat)

		tmp_drlon = np.gradient(tmp_rlon)
		tmp_drlat = np.gradient(tmp_rlat)

		tmp_dxi = np.sign(tmp_drlon[1])*const['a']*2*np.arcsin(np.abs(np.cos(tmp_rlat)*np.sin(tmp_drlon[1]*0.5)))
		tmp_dyi = np.sign(tmp_drlat[1])*const['a']*2*np.arcsin(np.abs(np.sin(tmp_drlat[1]*0.5)))
		tmp_di = (tmp_dxi**2+tmp_dyi**2)**0.5

		tmp_dxj = np.sign(tmp_drlon[0])*const['a']*2*np.arcsin(np.abs(np.cos(tmp_rlat)*np.sin(tmp_drlon[0]*0.5)))
		tmp_dyj = np.sign(tmp_drlat[0])*const['a']*2*np.arcsin(np.abs(np.sin(tmp_drlat[0]*0.5)))
		tmp_dj = (tmp_dxj**2+tmp_dyj**2)**0.5

		tmp_icurr = tmp_data['ucurr']*tmp_dxi/tmp_di+tmp_data['vcurr']*tmp_dyi/tmp_di
		tmp_jcurr = tmp_data['ucurr']*tmp_dxj/tmp_dj+tmp_data['vcurr']*tmp_dyj/tmp_dj

		###

		tmp_dtemp = np.gradient(tmp_data['temp'])	
		temp_ddepth = np.gradient(ncin_depth)							# [K]

		ncin_temp_ha = -(	tmp_icurr*(tmp_dtemp[3]/tmp_di)+ \
					tmp_jcurr*(tmp_dtemp[2]/tmp_dj)	\
				)

		ncin_temp_va = -(tmp_data['w'][:, 1:]+tmp_data['w'][:, :-1])/2*tmp_dtemp[1]/temp_ddepth[1]

		######

		np.savez('data/%s_%s/output_proc/ocean_tempadv.npz' % (tmp_tc, tmp_src, ),	hori=ncin_temp_ha, vert=ncin_temp_va,
												icurr=tmp_icurr, jcurr=tmp_jcurr,
												dtdi=(tmp_dtemp[3]/tmp_di), dtdj=(tmp_dtemp[2]/tmp_dj),
												dtdz=tmp_dtemp[1]/temp_ddepth[1],
												)
