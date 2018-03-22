# Introduction

This repository contains some classic or practical advection schemes for demonstrating them to first year graduates:

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

You can download it, and run each schemes by:

```
$ cd <repository directory>
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd upwind
$ ./upwind_adv_1d_case.exe ../../upwind/namelist
$ ncl ../../tools/plot_adv_1d.ncl scheme=\"upwind\"
```

# Software dependencies

- Fortran compiler (I use `gfortran`)
- NetCDF library (for outputting data)
- GSL library (for tridiagonal matrix solver)
- CMake (for generating Makefile)
- NCL (for plotting)

If you encounter difficulties to fulfill the requirement, you may use Docker image which we have created by:

```
$ docker pull dongli/iap-cgfd-adv-cases:0.0.1
$ docker run -it dongli/iap-cgfd-adv-cases:0.0.1
```

# Authors

- Li Dong <dongli@lasg.iap.ac.cn>
