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
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

######

plot_tc = ['HATO', 'NIDA', 'MANGKHUT', ]
plot_run = ['ctrl', 'r0rad', 'r0rad-to18', ]

######

def digitime(t):
	return np.datetime64(dt.datetime.strptime(t.decode('utf-8'), '%Y%m%d%H%M'))

def ms2knot(v):
	return float(v)*3.6/1.852

###

jttrk_dtype = np.dtype({'names':	('time',		'lon',	'lat',	'vm',	'cp', ),
			'formats':	('datetime64[m]',	'f4',	'f4',	'f4',	'f4', ),
			})

track_dtype = np.dtype([	('time', 'datetime64[m]'), ('minute', 'i2'), ('lon', 'f4'), ('lat', 'f4'), ('rawlon', 'f4'), ('rawlat', 'f4'),
				('cp', 'f4'), ('vm', 'f4'), ('vor850', 'f4'),
				('lndreg', 'i1'), ('lnddis', 'i1'),
				])

######

plot_mproj = plt_mcrs.PlateCarree()
plot_lon = np.arange(104, 128.1, 4)
plot_lat = np.arange(14, 26.1, 4)
plot_ntc = len(plot_tc)

###

plot_label = {	'ctrl':		'Control run',
		'r0rad':	'ROMS no radiation (NORAD)',
		'r0rad-to18':	'ROMS no radiation to 21/18Z (NORAD-t18)',
		'jtwc':		'JTWC best track',
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

cs_tidx = {	'HATO':		{	'ctrl':		43,
					'r0rad':	43,
					'r0rad-to18':	43,
					},
		'MANGKHUT':	{	'ctrl':		49,
					'r0rad':	49,
					},
		'NIDA':		{	'ctrl':		40,
					'r0rad':	40,
					},
		}
		
lnd_tidx = {	'HATO':		{	'ctrl':		57,
					'r0rad':	55,
					'r0rad-to18':	57,
					},
		'MANGKHUT':	{	'ctrl':		65,
					'r0rad':	65,
					},
		'NIDA':		{	'ctrl':		52,
					'r0rad':	52,
					},
		}
######

print('Opening NetCDF file: /home/clok/gwork/COAWST/making_scs_tc/wrf_grid_d01.nc')
ncdatain = Dataset('/home/clok/gwork/COAWST/making_scs_tc/wrf_grid_d01.nc')

nc_xlon = ncdatain.variables['XLONG'][0]
nc_xlat = ncdatain.variables['XLAT'][0]
nc_bathy = ncdatain.variables['BATHY'][0]

######

plt.rcParams['font.size'] = 12
fig = plt.figure(figsize=(21, 14))
ax = {0: {}, 1:{}}

###

for tmp_tcidx, tmp_tc in enumerate(plot_tc):
	tctrk = {}
	for tmp_src in plot_run:
		tmp_fntc = 'data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, )
		if not os.path.isfile(tmp_fntc) or not os.access(tmp_fntc, os.R_OK):
			continue

		tctrk[tmp_src] = np.loadtxt(tmp_fntc, dtype=track_dtype, ndmin=1, converters={7:ms2knot})
	tctrk['jtwc'] = np.loadtxt('data/%s_jtwc.txt' % (tmp_tc, ), dtype=jttrk_dtype, usecols=range(5), ndmin=1, converters={0:digitime})

	######

	ax[0][tmp_tcidx] = plt.subplot(2, 3, 1+tmp_tcidx, projection=plot_mproj)
	ax[0][tmp_tcidx].set_extent([plot_lon.min(), plot_lon.max(), plot_lat.min(), plot_lat.max()], crs=plot_mproj)

	for tmp_src, tmp_tctrk in tctrk.items():
		tmp_tctidx_0z = np.where(pd.DatetimeIndex(tmp_tctrk['time']).hour == 0)[0]
		img11 = ax[0][tmp_tcidx].plot(tmp_tctrk['lon'], tmp_tctrk['lat'], plot_ls[tmp_src], color=plot_colour[tmp_src], label='%s' % (plot_label[tmp_src], ))
		img12 = ax[0][tmp_tcidx].plot(tmp_tctrk['lon'][tmp_tctidx_0z], tmp_tctrk['lat'][tmp_tctidx_0z], plot_ms[tmp_src], color=plot_colour[tmp_src], ms=7)

	###

	img22 = ax[0][tmp_tcidx].contour(	nc_xlon, nc_xlat,
						np.abs(nc_bathy), [200],
						linestyles='dashed', colors='grey',
						)

	###

	ax[0][tmp_tcidx].set_title('(%s) Tracks of %s' % (chr(96+1+tmp_tcidx), tmp_tc.title()))

	ax[0][tmp_tcidx].set_xticks(plot_lon, crs=plot_mproj)
	ax[0][tmp_tcidx].xaxis.set_major_formatter(LongitudeFormatter())
	ax[0][tmp_tcidx].set_xlabel('Longitude')

	ax[0][tmp_tcidx].set_yticks(plot_lat, crs=plot_mproj)
	ax[0][tmp_tcidx].yaxis.set_major_formatter(LatitudeFormatter())
	if tmp_tcidx == 0:
		ax[0][0].set_ylabel('Latitude')

	ax[0][tmp_tcidx].grid()
	ax[0][tmp_tcidx].add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'), zorder=0)

	######

	ax[1][tmp_tcidx] = plt.subplot(2, 3, 4+tmp_tcidx)

	####

	for tmp_src, tmp_tctrk in tctrk.items():
		img11 = ax[1][tmp_tcidx].plot(tmp_tctrk['time'], tmp_tctrk['cp'], plot_ls[tmp_src], color=plot_colour[tmp_src], label='%s' % (plot_label[tmp_src], ))
		if 'jtwc' not in tmp_src:
			img12 = ax[1][tmp_tcidx].plot(	tmp_tctrk['time'][cs_tidx[tmp_tc][tmp_src]:lnd_tidx[tmp_tc][tmp_src]+1],
							tmp_tctrk['cp'][cs_tidx[tmp_tc][tmp_src]:lnd_tidx[tmp_tc][tmp_src]+1], plot_ms[tmp_src], color=plot_colour[tmp_src])

	###

	ax[1][tmp_tcidx].set_title('(%s) Intensity of %s' % (chr(96+4+tmp_tcidx), tmp_tc.title()))
	if tmp_tcidx == 0:
		ax[1][tmp_tcidx].legend()

	ax[1][tmp_tcidx].xaxis.set_major_locator(plt_dates.DayLocator())
	ax[1][tmp_tcidx].xaxis.set_major_formatter(plt_dates.DateFormatter('%d-%b'))
	ax[1][tmp_tcidx].xaxis.set_minor_locator(plt_dates.HourLocator(byhour=range(0, 24, 6)))
	ax[1][tmp_tcidx].set_xlim(tctrk['ctrl']['time'][0], tctrk['ctrl']['time'][-1])
	ax[1][tmp_tcidx].set_xlabel('UTC Time')

	if tmp_tcidx == 0:
		ax[1][0].set_ylabel('Central pressure / hPa')


	ax[1][tmp_tcidx].grid()

fig.subplots_adjust(hspace=0.1)
plt.savefig('ncee_%s_tctrk.png' % ('-'.join(plot_tc), ), bbox_inches='tight', dpi=300)
plt.close(fig)
