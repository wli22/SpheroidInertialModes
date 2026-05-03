# SpheroidInertialModes

This repository contains code for computation of inertial modes in spheres and spheroids.

We employ James Bremer's algorithm to compute the associated Legendre polynomials. The Fortran code was modified to exclude functions that compute the associated Legendre functions of the second kind in order to speed up the computation of inertial modes.

To use the code, you have to download five files, including alegendre_data.bin1, alegendre_data.bin2, alegendre_data.bin3, bessel_data.bin, and bessel_data16.bin from https://github.com/JamesCBremerJr/ALegendreEval
