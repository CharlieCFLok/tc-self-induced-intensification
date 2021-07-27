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
plot_run = ['ctrl', ]

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
plot_dep_max = 75

######

for tmp_tc in plot_tc:
	for tmp_src in plot_run:
		tmp_fntc = 'data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, )
		if not os.path.isfile(tmp_fntc) or not os.access(tmp_fntc, os.R_OK):
			continue
		tctrk = np.loadtxt(tmp_fntc, dtype=track_dtype, ndmin=1, converters={7:ms2knot})

		######

		tmp_data = {}
		
		for tmp_fnnc in os.listdir('data/%s_%s/output_nc/' % (tmp_tc, tmp_src, )):
			if 'ocean_his_d01' not in tmp_fnnc:
				continue

			print('Opening NetCDF file: data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))
			ncdatain = Dataset('data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))

			nc_xtime = ncdatain.variables['ocean_time']						# netcdf delta time

			nc_xlon = ncdatain.variables['lon_rho'][:]						# longitude on mass point [deg]
			nc_xlat = ncdatain.variables['lat_rho'][:]						# latitude on mass point [deg]

			nc_temp = ncdatain.variables['temp'][:]							# potential temperature on mass point [C]

			###

			tmp_tstripe = int(1800/(nc_xtime[1]-nc_xtime[0]))

			if 'xtime' in tmp_data:
				tmp_data['xtime'] = np.concatenate((tmp_data['xtime'], nc_xtime[tmp_tstripe::tmp_tstripe]), axis=None)
				tmp_data['temp'] = np.concatenate((tmp_data['temp'], nc_temp[tmp_tstripe::tmp_tstripe]), axis=0)
			else:
				tmp_data['xtime'] = nc_xtime[::tmp_tstripe]
				tmp_data['temp'] = nc_temp[::tmp_tstripe]

		###

		nc_time = pd.to_datetime(nc_xtime.units.split(' since ')[1])+pd.to_timedelta(tmp_data['xtime'], nc_xtime.units.split(' since ')[0])	# time
		nc_time_fix = tctrk['time'][0]-nc_time[0]-pd.to_timedelta(tctrk['minute'][0], 'minute')
		nc_time = nc_time+nc_time_fix

		nc_dtime = nc_time[1]-nc_time[0]							# time resolution

		nc_ntime, nc_ndepth, nc_nxlat, nc_nxlon = tmp_data['temp'].shape

		######

		npzdata_depth = np.load('data/%s_%s/output_proc/ocean_depth.npz' % (tmp_tc, tmp_src, ))
		nc_depth = npzdata_depth['depth'][:]

		######


		hovmo_d_lon = np.array([114.0, 115.0, 116.0, 115.5])
		hovmo_d_lat = np.array([21.0, 21.5, 22.0, 20.5])

		if tmp_tc == 'MANGKHUT':
			hovmo_d_lon = np.array([112.5, 113.5, 114.5, 114.0])
			hovmo_d_lat = np.array([20.5, 21.0, 21.5, 20.0])

		hovmo_d_temp = np.zeros(hovmo_d_lon.shape+(nc_ntime, nc_ndepth,))
		hovmo_d_xdepth = np.zeros_like(hovmo_d_temp)

		hovmo_d_xlon = np.zeros_like(hovmo_d_lon)
		hovmo_d_xlat = np.zeros_like(hovmo_d_lat)

		hovmo_d_xtime = np.broadcast_to(nc_time[:], (nc_ndepth, nc_ntime,)).T

		hovmo_d_cpa_tctidx = np.zeros(hovmo_d_lon.shape, dtype='int')

		######

		for tmp_hovmo_d_lonidx in range(hovmo_d_lon.shape[0]):
			tmp_dis = ((nc_xlon[:]-hovmo_d_lon[tmp_hovmo_d_lonidx])**2+(nc_xlat[:]-hovmo_d_lat[tmp_hovmo_d_lonidx])**2)**0.5

			hovmo_d_xlon[tmp_hovmo_d_lonidx] = nc_xlon[tmp_dis <= 0.25].mean()
			hovmo_d_xlat[tmp_hovmo_d_lonidx] = nc_xlat[tmp_dis <= 0.25].mean()

			for tmp_tidx in range(nc_time.shape[0]):
				for tmp_dlevidx in np.arange(nc_ndepth):
					tmp_temp = tmp_data['temp'][tmp_tidx, tmp_dlevidx].squeeze()
					tmp_xdepth = nc_depth[tmp_tidx, tmp_dlevidx].squeeze()

					hovmo_d_temp[tmp_hovmo_d_lonidx, tmp_tidx, tmp_dlevidx] = np.nanmean(tmp_temp[tmp_dis <= 0.25])
					hovmo_d_xdepth[tmp_hovmo_d_lonidx, tmp_tidx, tmp_dlevidx] = np.nanmean(tmp_xdepth[tmp_dis <= 0.25])

			######

			tmp_tcdis = ((tctrk['lon']-hovmo_d_lon[tmp_hovmo_d_lonidx])**2+(tctrk['lat']-hovmo_d_lat[tmp_hovmo_d_lonidx])**2)**0.5
			hovmo_d_cpa_tctidx[tmp_hovmo_d_lonidx] = np.where(tmp_tcdis == tmp_tcdis.min())[0]

			#####

		plt.rcParams['font.size'] = 12
		fig, axes = plt.subplots(2, 2, sharey=True, sharex=True, figsize=(18, 12))

		for tmp_hovmo_d_lonidx in range(hovmo_d_lon.shape[0]):
			img11 = axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].contourf(	hovmo_d_xtime, hovmo_d_xdepth[tmp_hovmo_d_lonidx],
									hovmo_d_temp[tmp_hovmo_d_lonidx], np.linspace(20, 31, num=23),
									extend='both', cmap='RdYlBu_r',
									)
			img12 = axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].contour(	hovmo_d_xtime, hovmo_d_xdepth[tmp_hovmo_d_lonidx],
									hovmo_d_temp[tmp_hovmo_d_lonidx], np.arange(18, 33, 2),
									colors='k',
									)

			img2 = axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].plot([tctrk['time'][hovmo_d_cpa_tctidx[tmp_hovmo_d_lonidx]], tctrk['time'][hovmo_d_cpa_tctidx[tmp_hovmo_d_lonidx]]], [0, -plot_dep_max], 'g--')

			###

			clabel = axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].clabel(img12, inline=1, fmt='%2d')

			axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].set_title('(%s) Location %s ' % (chr(96+1+tmp_hovmo_d_lonidx), chr(65+tmp_hovmo_d_lonidx), ))
			if int(tmp_hovmo_d_lonidx/2) == 1:
				axes[1, tmp_hovmo_d_lonidx%2].set_xlabel('UTC Time')

			axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].set_ylim(max(hovmo_d_xdepth.min(), -plot_dep_max), 0)
			if tmp_hovmo_d_lonidx%2 == 0:
				axes[int(tmp_hovmo_d_lonidx/2), 0].set_ylabel('Depth / m')

		cax = plt.axes([axes[-1, -1].get_position().x1+0.01, axes[-1, -1].get_position().y0, 0.01, axes[0, -1].get_position().y1-axes[-1, -1].get_position().y0])
		plt.colorbar(img11, cax=cax)

		axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].set_xlim([nc_time[0], nc_time[-1]])
		axes[int(tmp_hovmo_d_lonidx/2), tmp_hovmo_d_lonidx%2].xaxis.set_major_formatter(plt_dates.DateFormatter('%d/%HZ'))

		fig.subplots_adjust(wspace=0.1, hspace=0.1)
		plt.savefig('ncee_hover_%s_r-temp-%s.png' % (tmp_tc, tmp_src, ), bbox_inches='tight', dpi=300)
		plt.close(fig)
