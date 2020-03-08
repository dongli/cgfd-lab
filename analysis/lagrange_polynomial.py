#!/usr/bin/env python3

import numpy as np

def lag_poly_coef(x, x0):
	w = np.ones(len(x))
	for i in range(len(x)):
		for j in range(len(x)):
			if i != j:
				w[i] = w[i] * (x0 - x[j]) / (x[i] - x[j])
	return w

print(lag_poly_coef([-1.0, 0.0, 1.0, 2.0], 0.5))
print(lag_poly_coef([-1.0, 0.0, 1.0     ], 0.5))
print(lag_poly_coef([      0.0, 1.0, 2.0], 0.5))
