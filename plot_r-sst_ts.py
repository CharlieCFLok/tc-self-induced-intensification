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

plot_nature = False
if len(sys.argv) == 2:
	plot_nature = (sys.argv[1].lower() == 'nature')

######

plot_tc = ['HATO', 'NIDA', 'MANGKHUT' ]
plot_run = ['ctrl', 'r0rad', 'r0rad-to18', ]

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
plot_lon = np.arange(110, 120.1, 2)
plot_lat = np.arange(18, 24.1, 2)

###

plot_tsdepidx = -1

###

plot_label = {	'ctrl':		'Control',
		'r0rad':	'NORAD',
		'r0rad-to18':	'NORAD-t18',
		'jtwc':		'JTWC',
		}

plot_colour = {	'ctrl':		'blue',
		'r0rad':	'green',
		'r0rad-to18':	'purple',
		'jtwc':		'red',
		}
plot_ls = {	'ctrl':		'--',
		'r0rad':	':',
		'r0rad-to18':	'-.',
		'jtwc':		'-',
		}

plot_ms = {	'ctrl':		'o',
		'r0rad':	'd',
		'r0rad-to18':	'+',
		'jtwc':		's',
		}
######

for tmp_tc in plot_tc:
	plot_tslon = np.array([114.0, 115.0, 116.0, 115.5])
	plot_tslat = np.array([21.0, 21.5, 22.0, 20.5])

	if tmp_tc == 'MANGKHUT':
		plot_tslon = np.array([112.5, 113.5, 114.5, 114.0])
		plot_tslat = np.array([20.5, 21.0, 21.5, 20.0])

	plt.rcParams['font.size'] = 12
	fig =  plt.figure(figsize=(24, 12))
	gs = fig.add_gridspec(2, 3)

	######

	axes = {0: {}, 1: {}, }

	axes[0][0] = fig.add_subplot(gs[0, 0], projection=plot_mproj)
	axes[0][0].set_extent([plot_lon.min(), plot_lon.max(), plot_lat.min(), plot_lat.max()], crs=plot_mproj)

	######

	for tmp_src in plot_run:
		tmp_fntc = 'data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, )
		if not os.path.isfile(tmp_fntc) or not os.access(tmp_fntc, os.R_OK):
			continue
		tmp_tctrk = np.loadtxt(tmp_fntc, dtype=track_dtype, ndmin=1, converters={7:ms2knot})

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

			nc_bathy = ncdatain.variables['h'][:]

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
		nc_time_fix = tmp_tctrk['time'][0]-nc_time[0]-pd.to_timedelta(tmp_tctrk['minute'][0], 'minute')
		nc_time = nc_time+nc_time_fix

		nc_dtime = nc_time[1]-nc_time[0]							# time resolution

		nc_ntime, nc_ndepth, nc_nxlat, nc_nxlon = tmp_data['temp'].shape

		######

		npzdata_depth = np.load('data/%s_%s/output_proc/ocean_depth.npz' % (tmp_tc, tmp_src, ))
		nc_depth = npzdata_depth['depth'][:]

		######

		img01 = axes[0][0].plot(tmp_tctrk['lon'], tmp_tctrk['lat'], plot_ls[tmp_src], color=plot_colour[tmp_src], label='%s' % (plot_label[tmp_src], ))

		tmp_tctidx_0z = np.where(pd.DatetimeIndex(tmp_tctrk['time']).hour%12 == 0)[0]
		img011 = axes[0][0].plot(tmp_tctrk['lon'][tmp_tctidx_0z], tmp_tctrk['lat'][tmp_tctidx_0z], plot_ms[tmp_src], color=plot_colour[tmp_src], ms=8)

		if tmp_src == plot_run[0]:
			for tmp_tidx in tmp_tctidx_0z:
				if pd.to_datetime(tmp_tctrk['time'][tmp_tidx]).hour == 12:
					img012 = axes[0][0].text(	tmp_tctrk['lon'][tmp_tidx], tmp_tctrk['lat'][tmp_tidx]+0.2,
									pd.to_datetime(tmp_tctrk['time'][tmp_tidx]).strftime('%d/%HZ'),
									color=plot_colour[tmp_src], ha='center', va='bottom', clip_on=True)

			###

			img02 = axes[0][0].contour(	nc_xlon, nc_xlat,
							np.abs(nc_bathy), np.linspace(200, 1000, num=3),
							linestyles='dashed', colors='grey',
							)

		######

		for tmp_tsidx in range(plot_tslon.shape[0]):
			tmp_dis = ((nc_xlon[:]-plot_tslon[tmp_tsidx])**2+(nc_xlat[:]-plot_tslat[tmp_tsidx])**2)**0.5

			tmp_tslon = nc_xlon[tmp_dis <= 0.25].mean()
			tmp_tslat = nc_xlat[tmp_dis <= 0.25].mean()

			tmp_tsbathy = nc_bathy[tmp_dis <= 0.25].mean()

			tmp_tstemp = np.zeros(nc_ntime)
			tmp_tsdepth = np.zeros(nc_ntime)

			for tmp_tidx in range(nc_time.shape[0]):
				tmp_temp = tmp_data['temp'][tmp_tidx, plot_tsdepidx].squeeze()
				tmp_depth = nc_depth[tmp_tidx, plot_tsdepidx].squeeze()

				tmp_tstemp[tmp_tidx] = np.nanmean(tmp_temp[tmp_dis <= 0.25])
				tmp_tsdepth[tmp_tidx] = np.nanmean(tmp_depth[tmp_dis <= 0.25])

			######

			tmp_tcdis = ((tmp_tctrk['lon']-tmp_tslon)**2+(tmp_tctrk['lat']-tmp_tslat)**2)**0.5

			tmp_cpa_trktidx = np.where(tmp_tcdis == tmp_tcdis.min())[0][0]
			tmp_cpa_tidx = np.where(nc_time == pd.to_datetime(tmp_tctrk['time'][tmp_cpa_trktidx]))

			######

			if tmp_src == plot_run[0]:
				img031 = axes[0][0].plot(tmp_tslon, tmp_tslat, 'x', color='black')
				img032 = axes[0][0].text(tmp_tslon, tmp_tslat-0.2, chr(65+tmp_tsidx), color='black', ha='center', va='top')
				axes[0][0].set_title('(a) Tracks of %s' % (tmp_tc.title(), ))
			######

			tmp_axyidx = int((tmp_tsidx+1)/3)
			tmp_axxidx = (tmp_tsidx+1)%3+(plot_tslon.shape[0]==4)*(tmp_tsidx>1)

			if tmp_axxidx not in axes[tmp_axyidx].keys():
				axes[tmp_axyidx][tmp_axxidx] = fig.add_subplot(gs[tmp_axyidx, tmp_axxidx])

			img11 = axes[tmp_axyidx][tmp_axxidx].plot(nc_time, tmp_tstemp, plot_ls[tmp_src], color=plot_colour[tmp_src])
			img12 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_tctrk['time'][tmp_cpa_trktidx], tmp_tstemp[tmp_cpa_tidx], plot_ms[tmp_src], color=plot_colour[tmp_src], ms=8)

			if tmp_src == plot_run[0]:
				axes[tmp_axyidx][tmp_axxidx].set_title('(%s) Location %s (%5.1fE, %5.1fN, %3.0fm)' % (chr(96+2+tmp_tsidx), chr(64+1+tmp_tsidx), tmp_tslon, tmp_tslat, tmp_tsbathy, ))

				axes[tmp_axyidx][tmp_axxidx].set_xlim([nc_time[0], nc_time[-1]])
				axes[tmp_axyidx][tmp_axxidx].xaxis.set_major_formatter(plt_dates.DateFormatter('%d/%HZ'))
				if tmp_axyidx == 1:
					axes[1][tmp_axxidx].set_xlabel('UTC Time')

				if tmp_axxidx == 1:
					axes[tmp_axyidx][1].set_ylabel('SST / $^o$C')
				axes[tmp_axyidx][tmp_axxidx].yaxis.set_ticks(	np.arange(	np.floor(axes[tmp_axyidx][tmp_axxidx].get_ylim()[0]/0.5)*0.5,
												np.ceil(axes[tmp_axyidx][tmp_axxidx].get_ylim()[1]/0.5)*0.5+0.1,
												0.5,
												)
										)

	######

	axes[0][0].legend()

	axes[0][0].set_xticks(plot_lon, crs=plot_mproj)
	axes[0][0].set_xlabel('Longitude')

	axes[0][0].set_yticks(plot_lat, crs=plot_mproj)
	axes[0][0].set_ylabel('Latitude')

	axes[0][0].grid()
	axes[0][0].add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'))

	###

	plt.savefig('ncee_%s_r-sst_ts.png' % (tmp_tc, ), bbox_inches='tight', dpi=300)
	plt.close(fig)
