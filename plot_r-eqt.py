#!/home/clok/anaconda3/bin/python3

### system
import sys, os
### netcdf
from netCDF4 import Dataset
import wrf
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

from scipy import signal as sp_signal

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

plot_tslon = np.array([114.0, 115.0, 116.0, 115.5])
plot_tslat = np.array([21.0, 21.5, 22.0, 20.5])

plot_label = {  'ctrl':         'Control',
                'r0rad':        'NORAD',
                }

######

for tmp_tc in plot_tc:
	plt.rcParams['font.size'] = 12
	fig, axes = plt.subplots(4, 2, figsize=(16, 24))

	for tmp_src in plot_run:
		tmp_fntc = 'data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, )
		if not os.path.isfile(tmp_fntc) or not os.access(tmp_fntc, os.R_OK):
			continue
		tmp_tctrk = np.loadtxt(tmp_fntc, dtype=track_dtype, ndmin=1, converters={7:ms2knot})

		######

		tmp_dia = {}
		ncocndia = Dataset('data/%s_%s/output_nc/ocean_dia.nc' % (tmp_tc, tmp_src, ))

		tmp_dia['xtime'] = ncocndia.variables['ocean_time']								# netcdf delta time

		tmp_dia['xlon'] = ncocndia.variables['lon_rho'][:]								# longitude on mass point [deg]
		tmp_dia['xlat'] = ncocndia.variables['lat_rho'][:]								# latitude on mass point [deg]

		tmp_dia['tempdia_hadv'] = ncocndia.variables['temp_hadv'][:]
		tmp_dia['tempdia_vadv'] = ncocndia.variables['temp_vadv'][:]
		tmp_dia['tempdia_hdiff'] = ncocndia.variables['temp_hdiff'][:]
		tmp_dia['tempdia_vdiff'] = ncocndia.variables['temp_vdiff'][:]
		tmp_dia['tempdia_rate'] = ncocndia.variables['temp_rate'][:]

		###

		tmp_dia['time'] =	pd.to_datetime(tmp_dia['xtime'].units.split(' since ')[1])+ \
					pd.to_timedelta(tmp_dia['xtime'][:], tmp_dia['xtime'].units.split(' since ')[0])	# time
		tmp_timefix = tmp_tctrk['time'][0]-pd.to_timedelta(tmp_tctrk['minute'][0], 'minute')-tmp_dia['time'][0]
		tmp_dia['time'] = tmp_dia['time']+tmp_timefix

		ncdia_dtime = tmp_dia['time'][1]-tmp_dia['time'][0]								# time resolution
		ncdia_ntime, ncdia_ndepth, ncdia_nxlat, ncdia_nxlon = tmp_dia['tempdia_rate'].shape

		######

		tmp_data = {}

		for tmp_fnnc in os.listdir('data/%s_%s/output_nc/' % (tmp_tc, tmp_src, )):
			if 'ocean_his_d01' not in tmp_fnnc:
				continue

			print('Opening NetCDF file: data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))
			ncocn = Dataset('data/%s_%s/output_nc/%s' % (tmp_tc, tmp_src, tmp_fnnc, ))

			ncocn_xtime = ncocn.variables['ocean_time']										# netcdf delta time
			tmp_data['xlon'] = ncocn.variables['lon_rho'][:]									# longitude on mass point [deg]
			tmp_data['xlat'] = ncocn.variables['lat_rho'][:]									# latitude on mass point [deg]

			if 'xtime' in tmp_data:
				tmp_data['xtime'] = np.concatenate((tmp_data['xtime'], ncocn_xtime[1:]), axis=None)
				tmp_data['xdep'] = np.concatenate((tmp_data['xdep'], ncocn.variables['z_rho'][1:]), axis=0)

				tmp_data['temp'] = np.concatenate((tmp_data['temp'], ncocn.variables['temp'][1:]), axis=0)

				tmp_data['senhf'] = np.concatenate((tmp_data['senhf'], ncocn.variables['sensible'][1:]), axis=0)
				tmp_data['lathf'] = np.concatenate((tmp_data['lathf'], ncocn.variables['latent'][1:]), axis=0)
				tmp_data['lwrad'] = np.concatenate((tmp_data['lwrad'], ncocn.variables['lwrad'][1:]), axis=0)
				tmp_data['swrad'] = np.concatenate((tmp_data['swrad'], ncocn.variables['swrad'][1:]), axis=0)
			else:
				tmp_data['xtime'] = ncocn_xtime[:]
				tmp_data['xdep'] = ncocn.variables['z_rho'][:]									# depth on mass point [deg]

				tmp_data['temp'] = ncocn.variables['temp'][:]									# potential temperature on mass point [C]

				tmp_data['senhf'] = ncocn.variables['sensible'][:]								# surface sensible heat flux [W m**(-2)]
				tmp_data['lathf'] = ncocn.variables['latent'][:]								# surface latent heat flux [W m**(-2)]
				tmp_data['lwrad'] = ncocn.variables['lwrad'][:]									# net longwave radiation flux [W m**(-2)]
				tmp_data['swrad'] = ncocn.variables['swrad'][:]									# solar shortwave radiation flux [W m**(-2)]

		###

		tmp_data['time'] =	pd.to_datetime(ncocn_xtime.units.split(' since ')[1])+ \
					pd.to_timedelta(tmp_data['xtime'], ncocn_xtime.units.split(' since ')[0])		# time
		tmp_timefix = tmp_tctrk['time'][0]-pd.to_timedelta(tmp_tctrk['minute'][0], 'minute')-tmp_data['time'][0]
		tmp_data['time'] = tmp_data['time']+tmp_timefix


		ncocn_dtime = tmp_data['time'][1]-tmp_data['time'][0]								# time resolution
		ncocn_ntime, ncocn_ndepth, ncocn_nxlat, ncocn_nxlon = tmp_data['temp'].shape

		######

		npzdata_mld = np.load('data/%s_%s/output_proc/ocean_mld.npz' % (tmp_tc, tmp_src, ))
		ncocn_mld = npzdata_mld['mld'][:]
		ncocn_ismld = npzdata_mld['ismld'][:]

		######

		npzdata_tempadv = np.load('data/%s_%s/output_proc/ocean_tempadv.npz' % (tmp_tc, tmp_src, ))
		tmp_data['temp_ha'] = npzdata_tempadv['hori'][:]
		tmp_data['temp_va']  = npzdata_tempadv['vert'][:]

		tmp_sdata = {}

		for tmp_var in ['temp', 'temp_ha', 'temp_va', ]:
			tmp_sdata[tmp_var] = (tmp_data[tmp_var]*ncocn_ismld).sum(axis=1)/ncocn_ismld.sum(axis=1)

		for tmp_var in ['senhf', 'lathf', 'lwrad', 'swrad', ]:
			tmp_sdata[tmp_var] = -tmp_data[tmp_var]/ncocn_mld/1027/3850

		for tmp_var in ['tempdia_rate', 'tempdia_hadv', 'tempdia_vadv', 'tempdia_hdiff', 'tempdia_vdiff', ]:
			tmp_sdata[tmp_var] = (tmp_dia[tmp_var]*ncocn_ismld[1:-1:2]).sum(axis=1)/ncocn_ismld[1:-1:2].sum(axis=1)

		######

		for tmp_tsidx in range(plot_tslon.shape[0]):
			tmp_dis = ((tmp_data['xlon'][:]-plot_tslon[tmp_tsidx])**2+(tmp_data['xlat'][:]-plot_tslat[tmp_tsidx])**2)**0.5

			tmp_tslon = tmp_data['xlon'][tmp_dis <= 0.25].mean()
			tmp_tslat = tmp_data['xlat'][tmp_dis <= 0.25].mean()

			######

			tmp_ts = {}
	
			for tmp_var in [	'temp',
						'temp_ha', 'temp_va',
						'senhf', 'lathf', 'lwrad', 'swrad', 
						'tempdia_rate', 'tempdia_hadv', 'tempdia_vadv', 'tempdia_hdiff', 'tempdia_vdiff',
						]:
				tmp_ntime = ncdia_ntime if 'dia' in tmp_var else ncocn_ntime
				tmp_ts[tmp_var] = np.zeros(tmp_ntime)

				for tmp_tidx in range(tmp_ntime):
					tmp_tsdata = tmp_sdata[tmp_var][tmp_tidx]
					tmp_ts[tmp_var][tmp_tidx] = np.nanmean(tmp_tsdata[tmp_dis <= 0.25])

			###

			tmp_ts2 = {}
			for tmp_var in [	'temp_ha', 'temp_va',
						'senhf', 'lathf', 'lwrad', 'swrad',
						]:
				tmp_ts2[tmp_var] = tmp_ts[tmp_var][:-2:2]*0.25+tmp_ts[tmp_var][1:-1:2]*0.5+tmp_ts[tmp_var][2::2]*0.25

			tmp_ts2['qnet'] = tmp_ts2['senhf']+tmp_ts2['lathf']+tmp_ts2['lwrad']+tmp_ts2['swrad']
			tmp_ts2['dyn'] = tmp_ts2['temp_ha']+tmp_ts2['temp_va']

			tmp_ts2['temp_vd'] = tmp_ts['tempdia_vdiff']-tmp_ts2['qnet']
			tmp_ts2['diff'] = tmp_ts['tempdia_hdiff']+tmp_ts2['temp_vd']

			tmp_ts2['rhs'] = tmp_ts2['dyn']+tmp_ts2['diff']+tmp_ts2['qnet']

			tmp_ts2['time'] = tmp_data['time'][2::2]
	
			######

			tmp_tcdis = ((tmp_tctrk['lon']-tmp_tslon)**2+(tmp_tctrk['lat']-tmp_tslat)**2)**0.5

			tmp_cpa_trktidx = np.where(tmp_tcdis == tmp_tcdis.min())[0][0]

			######

			tmp_axyidx = tmp_tsidx
			tmp_axxidx = 1*(tmp_src == 'r0rad')

			###

			img11 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_ts2['time'], (tmp_ts['temp'][2::2]-tmp_ts['temp'][:-2:2])/3600*10**5, '-', color='red', lw=2, label='SST')

			img12 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_ts2['time'], (tmp_ts2['rhs'])*10**5, '--', color='blue', lw=2, label='R.H.S.')
			img13 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_ts2['time'], (tmp_ts2['qnet'])*10**5, '-.', color='green', label='Qnet')
			if tmp_src == 'ctrl':
				img14 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_ts2['time'], (tmp_ts2['swrad'])*10**5, ':', color='purple', label='S.W. Rad.')
			else:
				img14 = axes[tmp_axyidx][tmp_axxidx].plot(tmp_ts2['time'], (tmp_ts2['lathf'])*10**5, ':', color='tab:orange', label='Latent heat')

			###

			tmp_ymin, tmp_ymax = axes[tmp_axyidx][tmp_axxidx].get_ylim()
			img3 = axes[tmp_axyidx][tmp_axxidx].plot(	[tmp_tctrk['time'][tmp_cpa_trktidx], tmp_tctrk['time'][tmp_cpa_trktidx]],
									[tmp_ymin, tmp_ymax], '--', color='grey')

			###
			axes[tmp_axyidx][tmp_axxidx].set_title('(%s) Location %s (%s)' % (chr(96+1+tmp_tsidx*2+(tmp_src == 'r0rad')), chr(64+1+tmp_tsidx), plot_label[tmp_src], ))

			axes[tmp_axyidx][tmp_axxidx].set_xlim([tmp_data['time'][0], tmp_data['time'][-1]])
			axes[tmp_axyidx][tmp_axxidx].xaxis.set_major_formatter(plt_dates.DateFormatter('%d/%HZ'))
			if tmp_axyidx == 3:
				 axes[3][tmp_axxidx].set_xlabel('UTC Time')

			axes[tmp_axyidx][tmp_axxidx].set_ylim([tmp_ymin, tmp_ymax])
			if tmp_axxidx == 0:
				axes[tmp_axyidx][0].set_ylabel('SST Tendency / 10$^{-5}$ $^{\circ}$C s$^{-1}$')

	axes[tmp_axyidx][0].legend()
	axes[tmp_axyidx][1].legend()
	plt.savefig('ncee_%s_r-dTdt_ts.png' % (tmp_tc, ), bbox_inches='tight', dpi=300)
