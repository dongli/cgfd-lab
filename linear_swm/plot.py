#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import sys
from glob import glob

with PdfPages('linear_swm.pdf') as pdf:
	fig = plt.figure(figsize=(15, 5))
	ax1 = fig.add_subplot(1, 3, 1)
	ax2 = fig.add_subplot(1, 3, 2)
	ax3 = fig.add_subplot(1, 3, 3)
	for time_step in range(int(sys.argv[1]) + 1):
		file_name = f'linear_swm_a_grid.{str(time_step).zfill(3)}.nc'
		if os.path.isfile(file_name):
			f = Dataset(file_name, 'r')
			x = f.variables['x'][:]
			y = f.variables['y'][:]
			xi = f.variables['xi'][:]
			yi = f.variables['yi'][:]
			h = f.variables['h'][0,:]
			u = f.variables['u'][0,:]
			v = f.variables['v'][0,:]
			ax1.pcolormesh(x, y, h, cmap='bwr', vmin=-0.002, vmax=0.002, edgecolors='k')
			ax1.set_xticks([])
			ax1.set_yticks([])
			ax1.set_title('A-grid')
		file_name = f'linear_swm_b_grid.{str(time_step).zfill(3)}.nc'
		if os.path.isfile(file_name):
			f = Dataset(file_name, 'r')
			x = f.variables['x'][:]
			y = f.variables['y'][:]
			h = f.variables['h'][0,:]
			ax2.pcolormesh(x, y, h, cmap='bwr', vmin=-0.002, vmax=0.002, edgecolors='k')
			ax2.set_xticks([])
			ax2.set_yticks([])
			ax2.set_title('B-grid')
		file_name = f'linear_swm_c_grid.{str(time_step).zfill(3)}.nc'
		if os.path.isfile(file_name):
			f = Dataset(file_name, 'r')
			x = f.variables['x'][:]
			y = f.variables['y'][:]
			h = f.variables['h'][0,:]
			ax3.pcolormesh(x, y, h, cmap='bwr', vmin=-0.002, vmax=0.002, edgecolors='k')
			ax3.set_xticks([])
			ax3.set_yticks([])
			ax3.set_title('C-grid')
		pdf.savefig()
	plt.close()
