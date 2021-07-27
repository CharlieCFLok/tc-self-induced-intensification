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
import datetime as dt
import pandas as pd

### plotting
import matplotlib.pyplot as plt
import matplotlib.dates as plt_dates
import cartopy.crs as plt_mcrs
import cartopy.feature as plt_mfeature

######

plot_tc = ['HATO', ]
plot_run = ['ctrl', 'r0rad', ]

######

def ms2knot(v):
	return float(v)*3.6/1.852

###

track_dtype = np.dtype([	('time', 'datetime64[m]'), ('minute', 'i2'), ('lon', 'f4'), ('lat', 'f4'), ('rawlon', 'f4'), ('rawlat', 'f4'),
				('cp', 'f4'), ('vm', 'f4'), ('vor850', 'f4'),
				('lndreg', 'i1'), ('lnddis', 'i1'),
				])

######

plot_mproj = plt_mcrs.PlateCarree()
plot_lon = np.arange(105, 125.1, 4)
plot_lat = np.arange(14, 26.1, 4)
plot_ntc = len(plot_tc)

plot_tidx0 = -1
plot_tidx1 = 0

###

plot_label = {	'lathf':	'Latent heat flux',
		}

######

for tmp_tcidx, tmp_tc in enumerate(plot_tc):
	tctrk = {}
	hovmo_i_lathf = {}
	hovmo_i_tcxlon = {}
	hovmo_i_tccrx = {}

	for tmp_src in plot_run:
		tctrk[tmp_src] = np.loadtxt('data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, ), dtype=track_dtype, ndmin=1, converters={7:ms2knot})

		######

		print('Opening NetCDF file: data/%s_%s/output_nc/ocean_his_d01.nc' % (tmp_tc, tmp_src, ))
		ncdatain = Dataset('data/%s_%s/output_nc/ocean_his_d01.nc' % (tmp_tc, tmp_src, ))

		nc_xtime = ncdatain.variables['ocean_time']						# netcdf delta time
		nc_time = pd.to_datetime(nc_xtime.units.split(' since ')[1])+pd.to_timedelta(nc_xtime[:], nc_xtime.units.split(' since ')[0])	# time
		nc_time_fix = tctrk[tmp_src]['time'][0]-nc_time[0]-pd.to_timedelta(tctrk[tmp_src]['minute'][0], 'minute')
		nc_time = nc_time+nc_time_fix

		nc_dtime = nc_time[1]-nc_time[0]							# time resolution

		nc_xlon = ncdatain.variables['lon_rho'][:]						# longitude on mass point [deg]
		nc_xlat = ncdatain.variables['lat_rho'][:]						# latitude on mass point [deg]

		nc_temp = ncdatain.variables['temp'][:]							# potential temperature on mass point [C]

		nc_lathf = ncdatain.variables['latent'][:]						# surface latent heat flux [W m**(-2)]

		nc_ntime, nc_ndepth, nc_nxlat, nc_nxlon = nc_temp.shape

		######

		npzdata_mld = np.load('data/%s_%s/output_proc/ocean_mld.npz' % (tmp_tc, tmp_src, ))
		nc_mld = npzdata_mld['mld5'][:]

		######

		hovmo_i_xiidx = [115, 280]
		hovmo_i_xjidx = 260
		hovmo_i_dlevidx = -1

		hovmo_i_lathf[tmp_src] = nc_lathf[:, hovmo_i_xjidx-2:hovmo_i_xjidx+3, hovmo_i_xiidx[0]:hovmo_i_xiidx[1]+1].mean(axis=-2)

		hovmo_i_xlon = nc_xlon[hovmo_i_xjidx-2:hovmo_i_xjidx+3, hovmo_i_xiidx[0]:hovmo_i_xiidx[1]+1].mean(axis=-2)
		hovmo_i_xlat = nc_xlat[hovmo_i_xjidx-2:hovmo_i_xjidx+3, hovmo_i_xiidx[0]:hovmo_i_xiidx[1]+1].mean(axis=-2)

		###

		hovmo_i_tcxlon[tmp_src] = np.zeros_like(tctrk[tmp_src]['time'], dtype='f4')
		hovmo_i_tccrx[tmp_src] = np.zeros_like(tctrk[tmp_src]['time'], dtype='bool')

		for tmp_tctidx in range(hovmo_i_tcxlon[tmp_src].shape[0]):
			tmp_tcdis = ((tctrk[tmp_src]['lon'][tmp_tctidx]-hovmo_i_xlon)**2+(tctrk[tmp_src]['lat'][tmp_tctidx]-hovmo_i_xlat)**2)
			hovmo_i_tcxlon[tmp_src][tmp_tctidx] = hovmo_i_xlon[np.argmin(tmp_tcdis)]
			hovmo_i_tccrx[tmp_src][tmp_tctidx] = tmp_tcdis.min() < 0.1

	######

	plt.rcParams['font.size'] = 12
	fig = plt.figure(figsize=(7, 9))

	ax = plt.subplot(1, 1, 1)

	img11 = ax.contourf(	hovmo_i_xlon, nc_time[:],
				(hovmo_i_lathf['r0rad']-hovmo_i_lathf['ctrl']), np.linspace(-600, 600, num=21),
				extend='both', cmap='RdBu',
				)
	img12 = ax.contour(	hovmo_i_xlon, nc_time[:],
				(hovmo_i_lathf['r0rad']-hovmo_i_lathf['ctrl']), np.arange(-600, 601, 240),
				colors='k',
				)

	img21 = ax.plot(hovmo_i_tcxlon['r0rad'], tctrk['r0rad']['time'], 'g-')
	img22 = ax.plot(hovmo_i_tcxlon['r0rad'][hovmo_i_tccrx['r0rad']], tctrk['r0rad']['time'][hovmo_i_tccrx['r0rad']], 'gs')

	###

	plt.colorbar(img11)
	clabel = ax.clabel(img12, inline=1, fmt='%2d')

	ax.set_xlim(hovmo_i_xlon[[0, -1]])
	ax.set_xlabel('Longitude')

	ax.set_ylim([nc_time[plot_tidx0], nc_time[plot_tidx1]])
	ax.yaxis.set_major_formatter(plt_dates.DateFormatter('%d/%HZ'))
	ax.set_ylabel('UTC Time')

	###

	plt.savefig('ncee_%s_r-lathf_ctrl-r0rad_hovmo-i.png' % (tmp_tc, ), bbox_inches='tight', dpi=300)
	plt.close(fig)
