# OpenQSWalk
Quantum stochastic walk simulation in open systems. 

Currently this repository consists of the program "ExpoAction" which calculates the action of the matrix exponential on a vector. The algorithm used is the most recent version of the Pade-approximate method developed by Al-Mohy et. al. in their 2009 paper "A New Scaling and Squaring Algorithm for the Matrix Exponential" and then updated in 2015. This algorithm is also used by the "SciPy" library.

The program is designed to work with compressed sparse row (CSR) matrices (the Harwell-Boeing format is supported for import and export) has been parallelised using OpenMp. This mostly occurs in Nullarbor.f90. “Nullarbor” is a library supporting CSR matrix operations import and export. Its development was necessitated by the requirements of the assignment for which “ExpoAction” was coded. 

High on the To-Do list is the replacement of Nullarbor with a (presumably) more optimised library (i.e. Intel MKL Sparse BLAS). 
