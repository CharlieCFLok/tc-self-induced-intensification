#!/home/clok/anaconda3/bin/python3

### system
#import sys, os
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

from matplotlib.colors import ListedColormap

######

plot_tc = 'HATO'
plot_run = 'ctrl'

plot_lon = np.arange(95, 145.1, 5)
plot_lat = np.arange(0, 40.1, 5)

######

datadir = 'data/%s_%s/output_nc/' % (plot_tc, plot_run, )

fnncin = {	'roms':	'ocean_his_d01.nc',
		'wrf1':	'wrfout_d01.nc',
		'wrf2':	'wrfout_d02.nc',
		}

###

nc_varname = {	'lon':	{	'roms':	'lon_rho',
				'wrf1':	'XLONG',
				'wrf2':	'XLONG',
				},
		'lat':	{	'roms':	'lat_rho',
				'wrf1':	'XLAT',
				'wrf2':	'XLAT',
				},
		}

###

plot_mproj = plt_mcrs.PlateCarree()

plot_colour = {	'wrf1':	'red',
		'wrf2':	'green',
		'roms':	'blue',
		}

######

plt.rcParams['font.size'] = 12
fig = plt.figure(figsize=(9, 7))

ax = plt.axes(projection=plot_mproj)
ax.set_extent([92, 145, 0, 40], crs=plot_mproj)

###

for tmp_model, tmp_fnncin in fnncin.items():
	print('Opening NetCDF file: '+datadir+'/'+tmp_fnncin)
	ncdatain = Dataset(datadir+'/'+tmp_fnncin)

	nc_xlon = ncdatain.variables[nc_varname['lon'][tmp_model]]	# longitude on mass point [deg]
	nc_xlat = ncdatain.variables[nc_varname['lat'][tmp_model]]	# latitude on mass point [deg]

	if len(nc_xlon[:].shape) == 2:
		plot_xlon = nc_xlon[:]
		plot_xlat = nc_xlat[:]
	elif len(nc_xlon[:].shape) == 3:
		plot_xlon = nc_xlon[0]
		plot_xlat = nc_xlat[0]

	###

	plot_data = np.zeros_like(plot_xlon)
	plot_data[2:-2, 2:-2] = np.nan


	img11 = ax.contourf(	plot_xlon, plot_xlat,
				plot_data,
				colors=plot_colour[tmp_model], zorder=10,
				)

###

print('Opening NetCDF file: data/wp.nc')
ncdatain = Dataset('data/wp.nc')

nc_lon = ncdatain.variables['x'][:]	# longitude on mass point [deg]
nc_lat = ncdatain.variables['y'][:]
nc_bathy = ncdatain.variables['z'][:]

plot_cmap = ListedColormap(('skyblue', 'greenyellow', ))
img12 = ax.contourf(	nc_lon, nc_lat,
			nc_bathy, [-1000, -200, 10, ],
			cmap=plot_cmap, zorder=1,
			)

img21 = ax.plot([111.18, 120.06], [19.77, 22.52], 'k--')

######

ax.set_xticks(plot_lon, crs=plot_mproj)
ax.set_xlabel('Longitude')

ax.set_yticks(plot_lat, crs=plot_mproj)
ax.set_ylabel('Latitude')

ax.grid()
ax.add_feature(plt_mfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='lightgrey'), zorder=5)

plt.savefig('ncee_domain_grid.png' , bbox_inches='tight', dpi=300)
plt.close(fig)
