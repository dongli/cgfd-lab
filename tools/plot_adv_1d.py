#!/usr/bin/env python3

import argparse
from glob import glob
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import sys

parser = argparse.ArgumentParser(description='Plot 1d advection test cases.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(      dest='file_patterns', help='NetCDF esult file pattern', nargs='+')
parser.add_argument('-s', dest='time_step_stride', help='Time step stride', default=1, type=int)
parser.add_argument('-e', dest='end_time_step', help='End time step', type=int)
args = parser.parse_args()

colors = ['black', 'red', 'blue', 'green', 'purple', 'magenta']

file_paths = []
for file_pattern in args.file_patterns:
	file_paths.append(sorted(glob(file_pattern)))

pdf_file_name = '_vs_'.join([os.path.basename(file_pattern).split('.')[0] for file_pattern in args.file_patterns])

with PdfPages(pdf_file_name + '.pdf') as pdf:
	for k in range(0, len(file_paths[0]), args.time_step_stride):
		fig = plt.figure(figsize=(8, 8))
		ax = fig.add_subplot(1, 1, 1)
		for j in range(len(file_paths)):
			f = Dataset(file_paths[j][k], 'r')
			ax.plot(f.variables['x'][:], f.variables['rho'][0,:], color=colors[j], linewidth=3, label=f.scheme.upper())
			f.close()
		ax.legend()
		ax.set_xlim(0, 1)
		ax.set_ylim(-0.1, 1.1)
		ax.tick_params(labelsize=18)
		pdf.savefig()
		plt.close()
