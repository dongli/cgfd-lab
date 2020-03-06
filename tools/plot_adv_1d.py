#!/usr/bin/env python3

import argparse
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import MFDataset
import os
import sys

parser = argparse.ArgumentParser(description='Plot 1d advection test cases.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', dest='file_pattern', help='NetCDF esult file pattern', required=True)
parser.add_argument('-s', dest='time_step_stride', help='Time step stride', default=1, type=int)
parser.add_argument('-e', dest='end_time_step', help='End time step', type=int)
args = parser.parse_args()

f = MFDataset(args.file_pattern, 'r')

with PdfPages(os.path.basename(args.file_pattern).split('.')[0] + '.pdf') as pdf:
	for k in range(0, f.dimensions['time'].dimtotlen, args.time_step_stride):
		fig = plt.figure(figsize=(15, 8))
		ax = fig.add_subplot(1, 1, 1)
		ax.plot(f.variables['x'][:], f.variables['rho'][k,:], color='black', linewidth=3)
		ax.set_title(f.scheme, fontsize=18)
		ax.set_xlim(0, 1)
		ax.set_ylim(-0.1, 1.1)
		ax.tick_params(labelsize=18)
		pdf.savefig()
		plt.close()
