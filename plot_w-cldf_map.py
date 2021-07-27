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

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

######

plot_tc = ['HATO', 'MANGKHUT', ]
plot_run = ['ctrl', ]

plot_domain = 1

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

######

plt.rcParams['font.size'] = 12
fig = plt.figure(figsize=(16, 12))

for tmp_src in plot_run:
	for tmp_tcidx, tmp_tc in enumerate(plot_tc):
		tctrk = np.loadtxt('data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, tmp_src, ), dtype=track_dtype, ndmin=1, converters={7:ms2knot})

		######

		tmp_tidx = 30
		if plot_domain == 1:
			tmp_tidx = int(tmp_tidx/3)
		tmp_varidx = 0
		ptmp_var = 0

		###

		print('Opening NetCDF file: data/%s_%s/output_nc/wrfout_d%02d.nc' % (tmp_tc, tmp_src, plot_domain, ))
		ncdatain = Dataset('data/%s_%s/output_nc/wrfout_d%02d.nc' % (tmp_tc, tmp_src, plot_domain, ))

		nc_time = wrf.extract_times(ncdatain, tmp_tidx)					# time

		nc_cldf = wrf.getvar(ncdatain, 'cloudfrac', tmp_tidx)				
		nc_w = wrf.getvar(ncdatain, 'wa', tmp_tidx)

		nc_p = wrf.getvar(ncdatain, 'pressure', tmp_tidx)
		nc_w500 = wrf.interplevel(nc_w, nc_p, (500, 600, 700, 800, 900))

		nc_xlat, nc_xlon = wrf.latlon_coords(nc_cldf)

		######

		tmp_tctidx = np.where(tctrk['time'] == nc_time)[0]

		###

		ax = plt.subplot(plot_ntc, 2, 1+tmp_tcidx*2, projection=plot_mproj)

		tmp_blues = ListedColormap(cm.get_cmap('Blues', 512)(np.linspace(0, 0.75, 256)))
		img11 = ax.contourf(	nc_xlon, nc_xlat,
					nc_cldf[0], np.linspace(0.5, 1, num=21),
					extend='min', cmap=tmp_blues,
					)

		img12 = ax.contourf(	nc_xlon, nc_xlat,
					nc_cldf[1], [0.8, 1],
					colors='none', hatches=['---', ],
					)
		img12.collections[0].set_edgecolor('navy')
		img12.collections[0].set_linewidth(0)

		img13 = ax.contourf(	nc_xlon, nc_xlat,
					nc_cldf[2], [0.8, 1],
					colors='none', hatches=['|||', ],
					)
		img13.collections[0].set_edgecolor('red')
		img13.collections[0].set_linewidth(0)

		img21 = ax.plot(tctrk['lon'], tctrk['lat'], '-', lw=3, color='white')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx], tctrk['lat'][tmp_tctidx], 's', ms=8, color='white')

		circlex = np.cos(np.radians(np.linspace(0, 360, 361)))
		circley = np.sin(np.radians(np.linspace(0, 360, 361)))
		for tmp_ridx in [2, 4, 6]:
			img23 = ax.plot(tctrk['lon'][tmp_tctidx]+circlex*tmp_ridx, tctrk['lat'][tmp_tctidx]+circley*tmp_ridx, '-', lw=3, color='white')
			img24 = ax.plot(tctrk['lon'][tmp_tctidx]+circlex*tmp_ridx, tctrk['lat'][tmp_tctidx]+circley*tmp_ridx, '--', color='green')

		img21 = ax.plot(tctrk['lon'], tctrk['lat'], '-', color='green')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx], tctrk['lat'][tmp_tctidx], 's', color='green')

		if tmp_tcidx == plot_ntc-1:
			cax = plt.axes([ax.get_position().x0, ax.get_position().y0-0.06, ax.get_position().width, 0.01])
			plt.colorbar(img11, cax=cax, orientation='horizontal')

		ax.set_title('(%s) Cloud cover of %s' % (chr(96+1+tmp_tcidx*2), tmp_tc.title()))

		ax.set_xlim(plot_lon[[0, -1]])
		ax.set_xticks(plot_lon, crs=plot_mproj)
		if tmp_tcidx == plot_ntc-1:
			ax.set_xlabel('Longitude')

		ax.set_ylim(plot_lat[[0, -1]])
		ax.set_yticks(plot_lat, crs=plot_mproj)
		ax.set_ylabel('Latitude')

		ax.grid()
		ax.coastlines('50m')

		###

		ax = plt.subplot(plot_ntc, 2, 2+tmp_tcidx*2, projection=plot_mproj)
		img11 = ax.contourf(	nc_xlon, nc_xlat,
					np.nanmean(nc_w500, axis=0), np.linspace(-0.5, 0.5, num=21),
					extend='both', cmap='bwr',
					)

		img21 = ax.plot(tctrk['lon'], tctrk['lat'], '-', color='green')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx], tctrk['lat'][tmp_tctidx], 's', color='green')

		if tmp_tcidx == plot_ntc-1:
			cax = plt.axes([ax.get_position().x0, ax.get_position().y0-0.06, ax.get_position().width, 0.01])
			plt.colorbar(img11, cax=cax, orientation='horizontal')

		ax.set_title('(%s) 500-900 hPa ertical motion of %s' % (chr(96+2+tmp_tcidx*2), tmp_tc.title()))

		ax.set_xlim(plot_lon[[0, -1]])
		ax.set_xticks(plot_lon, crs=plot_mproj)
		if tmp_tcidx == plot_ntc-1:
			ax.set_xlabel('Longitude')

		ax.set_ylim(plot_lat[[0, -1]])
		ax.set_yticks(plot_lat, crs=plot_mproj)

		ax.grid()
		ax.coastlines('50m')

		

	fig.subplots_adjust(wspace=0.1, hspace=0.05)
	plt.savefig('ncee_%s_w-cldf_%s_t%03d.png' % ('-'.join(plot_tc), tmp_src, tmp_tidx, ), bbox_inches='tight', dpi=300)
	plt.close(fig)
