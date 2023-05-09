# Constrained Optimization

This folder contains the routines that implement various strategies for 
approximating the pth root of a stochastic matrix using constrained 
optimization:
- `approximatepower.m` implements the interior-point method optimization
for the objective function `f(x) = 0.5*||x^p - a||_F^2`;
- `approximateroot.m` implements the interior-point method optimization
for the objective function `f(x) = 0.5*||x - a^{1/p}||_F^2`;
constrained on the set of stochastic matrices.


