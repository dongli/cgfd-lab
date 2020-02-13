# Introduction

This repository contains some classic or practical advection schemes for demonstrating them to the first year graduates:

- FTCS (Forward-in-time, centered-in-space)
- Upwind
- Beam-Warming
- Lax-Wendroff
- Fromm
- Leap-frog
- Crank-Nicolson
- TSPAS
- MPDATA
- WENO
- FFSL
- Semi-Lagrangian

You can download it, and run each schemes by:

```
$ cd <repository directory>
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd advection/upwind
$ ./upwind_adv_1d_case.exe
$ ncl ../../../tools/plot_adv_1d.ncl scheme=\"upwind\"
```

# Software dependencies

- Fortran compiler (I use `gfortran`)
- NetCDF library (for outputting data)
- GSL library (optional for tridiagonal matrix solver)
- CMake (for generating Makefile)
- NCL (for plotting)

# TODO list

Next I would like to do the following things:

- Add 2D upwind, Lax-Wendroff scheme cases
- Extend TSPAS to 2D and include cross terms

# Authors

- Li Dong <dongli@lasg.iap.ac.cn>
