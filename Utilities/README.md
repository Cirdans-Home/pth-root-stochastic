# Auxiliary Routines

This directory contains auxiliary code which is used either in optimization routines, to calculate auxiliary quantities or to generate test problems.

- `pvgth.m` computes the steady state vector x of the stochastic matrix P through LU factorization of I-P, with diagonal adjustment.
- `perfprof.m` produces a peformance profile for the data in the M-by-N matrix A, where A(i,j) > 0 measures the performance of the j'th solver on the i'th problem