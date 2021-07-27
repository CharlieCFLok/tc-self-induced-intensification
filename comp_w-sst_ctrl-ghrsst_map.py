#!/home/clok/anaconda3/bin/python3

### system
#import sys, os
### netcdf
from netCDF4 import Dataset
import wrf
### numerical
import numpy as np
from scipy import interpolate as sp_intpl
### calendar
import datetime as dt
import pandas as pd

### plotting
import matplotlib.pyplot as plt
#import matplotlib.dates as plt_dates
import cartopy.crs as plt_mcrs
import cartopy.feature as plt_mfeature

######

plot_tc = ['HATO', ]
plot_src = 'ctrl'
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
plot_lon = np.arange(105, 125.1, 2)
plot_lat = np.arange(14, 26.1, 2)

######

for tmp_tc in plot_tc:
	tctrk = np.loadtxt('data/%s_%s/output_proc/tctrk.wrfout_d02.nc.tc01.txt' % (tmp_tc, plot_src, ), dtype=track_dtype, ndmin=1, converters={7:ms2knot})

	######

	print('Opening NetCDF file: data/%s_%s/output_nc/wrfout_d%02d.nc' % (tmp_tc, plot_src, plot_domain, ))
	ncwrfdatain = Dataset('data/%s_%s/output_nc/wrfout_d%02d.nc' % (tmp_tc, plot_src, plot_domain, ))

	ncwrf_time = wrf.extract_times(ncwrfdatain, wrf.ALL_TIMES)					# time
	ncwrf_xtime = wrf.extract_times(ncwrfdatain, wrf.ALL_TIMES, do_xtime=True)			# delta time from model initalisation [minute]
	ncwrf_dtime = ncwrf_time[1]-ncwrf_time[0]							# time resolution

	ncwrf_land = wrf.getvar(ncwrfdatain, 'XLAND', wrf.ALL_TIMES)					# land mask [1 = land, 2 = water]
	ncwrf_sst = wrf.getvar(ncwrfdatain, 'SST', wrf.ALL_TIMES)					# sea surface temperature [K]

	ncwrf_xlat, ncwrf_xlon = wrf.latlon_coords(ncwrf_land)						# latitude, longitude [deg]

	######

	tmp_wsst = wrf.to_np(ncwrf_sst)-273.15
	tmp_wsst[np.where(ncwrf_land == 1)] = np.nan

	tmp_wntime_day = (np.timedelta64(24, 'h')/ncwrf_dtime).astype('int')

	######

	plt.rcParams['font.size'] = 12
	fig = plt.figure(figsize=(27, 14))

	######

	for tmp_tidx in [1, 2]:
		tmp_datestr = pd.DatetimeIndex(ncwrf_time).strftime('%Y%m%d')[tmp_tidx*tmp_wntime_day]
		tmp_twsst = tmp_wsst[tmp_tidx*tmp_wntime_day:(tmp_tidx+1)*tmp_wntime_day].mean(axis=0)

		######

		print('Opening NetCDF file: data/ghrsst_abom/regrid_upk_sst_%s.nc' % (tmp_datestr, ))
		ncghrdatain = Dataset('data/ghrsst_abom/regrid_upk_sst_%s.nc' % (tmp_datestr, ))
		ncghrdatain.set_auto_mask(False)

		ncghr_lon = ncghrdatain.variables['lon']							# longitude on mass point [deg]
		ncghr_lat = ncghrdatain.variables['lat']							# latitude on mass point [deg]
		ncghr_xlon, ncghr_xlat = np.meshgrid(ncghr_lon, ncghr_lat)

		ncghr_sst = ncghrdatain.variables['analysed_sst'][:].squeeze()-273.15				# sea surface temperature [K]
		ncghr_mask = ncghrdatain.variables['mask'][:].squeeze()
		ncghr_sst[ncghr_mask > 1] = np.nan

		######

		ncwrf_lonlat_pnt = np.array((wrf.to_np(ncwrf_xlon).flatten(), wrf.to_np(ncwrf_xlat).flatten())).T

		tmp_twsst_ghr = sp_intpl.griddata(ncwrf_lonlat_pnt, tmp_twsst.flatten(), (ncghr_xlon, ncghr_xlat), method='linear')
		tmp_twsst_err = tmp_twsst_ghr-ncghr_sst
		print('ERROR RMS:  %20.5f' % (np.nanmean(tmp_twsst_err**2)**0.5))
		print('ERROR BIAS: %20.5f' % (np.nanmean(tmp_twsst_err)))

		######

		plot_ncol = 3 

		###

		ax = plt.subplot(2, plot_ncol, 1+(tmp_tidx-1)*plot_ncol, projection=plot_mproj)
		ax.set_extent([plot_lon.min(), plot_lon.max(), plot_lat.min(), plot_lat.max()], crs=plot_mproj)

		img11 = ax.contourf(	ncwrf_xlon, ncwrf_xlat,
					tmp_twsst, np.linspace(24, 32, num=17),
					extend='both', cmap='RdYlBu_r',
					)
		img12 = ax.contour(	ncwrf_xlon, ncwrf_xlat,
					tmp_twsst, np.arange(24, 33, 2),
					colors='k',
					)

		tmp_tctidx_day = np.where(pd.DatetimeIndex(tctrk['time']).day == pd.DatetimeIndex(ncwrf_time).day[tmp_tidx*tmp_wntime_day])[0]
		img21 = ax.plot(tctrk['lon'], tctrk['lat'], 'g-')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx_day], tctrk['lat'][tmp_tctidx_day], 'gs')

		###

		if tmp_tidx == 2:
			cax = plt.axes([ax.get_position().x0, ax.get_position().y0-0.08, ax.get_position().width, 0.02])
			plt.colorbar(img11, cax=cax, orientation='horizontal')
		clabel = ax.clabel(img12, inline=1, fmt='%2d')

		ax.set_title('(%s) Control run (%s)' % (chr(97+(tmp_tidx-1)*plot_ncol), pd.DatetimeIndex(ncwrf_time).strftime('%d-%b-%Y')[tmp_tidx*tmp_wntime_day], ))

		ax.set_xticks(plot_lon, crs=plot_mproj)
		if tmp_tidx == 2:
			ax.set_xlabel('Longitude')

		ax.set_yticks(plot_lat, crs=plot_mproj)
		ax.set_ylabel('Latitude')

		ax.grid()
		ax.add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'), zorder=0)

		######

		ax = plt.subplot(2, plot_ncol, 2+(tmp_tidx-1)*plot_ncol, projection=plot_mproj)
		ax.set_extent([plot_lon.min(), plot_lon.max(), plot_lat.min(), plot_lat.max()], crs=plot_mproj)

		img11 = ax.contourf(	ncghr_xlon, ncghr_xlat,
					ncghr_sst[:], np.linspace(24, 32, num=17),
					extend='both', cmap='RdYlBu_r',
					)
		img12 = ax.contour(	ncghr_xlon, ncghr_xlat,
					ncghr_sst[:], np.arange(24, 33, 2),
					colors='k',
					)

		tmp_tctidx_day = np.where(pd.DatetimeIndex(tctrk['time']).day == pd.DatetimeIndex(ncwrf_time).day[tmp_tidx*tmp_wntime_day])[0]
		img21 = ax.plot(tctrk['lon'], tctrk['lat'], 'g-')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx_day], tctrk['lat'][tmp_tctidx_day], 'gs')

		###

		if tmp_tidx == 2:
			cax = plt.axes([ax.get_position().x0, ax.get_position().y0-0.08, ax.get_position().width, 0.02])
			plt.colorbar(img11, cax=cax, orientation='horizontal')
		clabel = ax.clabel(img12, inline=1, fmt='%2d')

		ax.set_title('(%s) GHRSST (%s)' % (chr(97+1+(tmp_tidx-1)*plot_ncol), pd.DatetimeIndex(ncwrf_time).strftime('%d-%b-%Y')[tmp_tidx*tmp_wntime_day], ))

		ax.set_xticks(plot_lon, crs=plot_mproj)
		if tmp_tidx == 2:
			ax.set_xlabel('Longitude')

		ax.set_yticks(plot_lat, crs=plot_mproj)

		ax.grid()
		ax.add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'), zorder=0)

		######

		ax = plt.subplot(2, plot_ncol, 3+(tmp_tidx-1)*plot_ncol, projection=plot_mproj)
		ax.set_extent([plot_lon.min(), plot_lon.max(), plot_lat.min(), plot_lat.max()], crs=plot_mproj)

		img11 = ax.contourf(	ncghr_xlon, ncghr_xlat,
					tmp_twsst_err[:], np.linspace(-4, 4, num=17),
					extend='both', cmap='RdBu_r',
					)
		img12 = ax.contour(	ncghr_xlon, ncghr_xlat,
					tmp_twsst_err[:], np.arange(-3, 4, 2),
					colors='k',
					)

		tmp_tctidx_day = np.where(pd.DatetimeIndex(tctrk['time']).day == pd.DatetimeIndex(ncwrf_time).day[tmp_tidx*tmp_wntime_day])[0]
		img21 = ax.plot(tctrk['lon'], tctrk['lat'], 'g-')
		img22 = ax.plot(tctrk['lon'][tmp_tctidx_day], tctrk['lat'][tmp_tctidx_day], 'gs')

		###

		if tmp_tidx == 2:
			cax = plt.axes([ax.get_position().x0, ax.get_position().y0-0.08, ax.get_position().width, 0.02])
			plt.colorbar(img11, cax=cax, orientation='horizontal')
		clabel = ax.clabel(img12, inline=1, fmt='%2d')

		ax.set_title('(%s) %s-%s' % (chr(97+2+(tmp_tidx-1)*plot_ncol), chr(97+(tmp_tidx-1)*plot_ncol), chr(97+1+(tmp_tidx-1)*plot_ncol), ))

		ax.set_xticks(plot_lon, crs=plot_mproj)
		if tmp_tidx == 2:
			ax.set_xlabel('Longitude')

		ax.set_yticks(plot_lat, crs=plot_mproj)

		ax.grid()
		ax.add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'), zorder=0)

		###

	fig.subplots_adjust(wspace=0.1, hspace=0.05)

	plt.savefig('ncee_%s_w-sst_%s-ghr.png' % (tmp_tc, plot_src, ), bbox_inches='tight', dpi=300)
	plt.close(fig)
