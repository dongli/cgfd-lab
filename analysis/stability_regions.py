#!/usr/bin/env python3

import matplotlib
# matplotlib.use('agg')
# from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

def draw_stability_region(x, y, R, color):
	zlevel = np.abs(R)
	plt.contour(x, y, zlevel, (1, 1.00001), colors=(color))

x, y = np.meshgrid(np.linspace(-5, 5, 100), np.linspace(-5, 5, 100))
z = x + 1j * y

# 1st-order Euler forward
R = 1 + z
draw_stability_region(x, y, R, 'red')

# 2nd-order predictor-corrector
R = 1 + z + z**2 / 2 + z**3 / 4
draw_stability_region(x, y, R, 'blue')

# 4th-order Runge-Kutta
R = 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24
draw_stability_region(x, y, R, 'green')

R = 1 + z + z**2 / 2 + z**3 / 6
draw_stability_region(x, y, R, 'black')

plt.grid(b=True, which='major', color='gray', linestyle='--')
plt.show()
