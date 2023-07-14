%% Optimization with stochastic matrices with assigned pi vector
% Performance tests for the solution of the approximation problem on the
% set of row-stochastic matrices with assigned stationary vector. Comparison
% between the usage of the Riemannian optimization algorithms and the
% standard constrained optimization one.

clear all; clc; close all hidden;

try manopt_version;
    fprintf("Using MANOPT version %d.%d.%d\n\n",manopt_version);
catch
    fprintf("Loading MANOPT\n");
    here = pwd;
    cd ../manopt
    addpath('auxiliaries/');
    importmanopt;
    cd(here);
end
addpath('../Utilities/'); % code for performance profile plots
addpath('../Matrices/'); % code for the generation of test matrices
addpath('../ConstrainedOptimization/'); % code for constrained optimization
addpath('../Manifold'); % code for the new manifold
%% Generation of the test problems
seed = 42;  % Random number generator seed use it to enforce reproducibility
p = 2;      % Square-root, change it to look for other roots
%% Running optimization tests
fprintf("Riemannian trust region      & Riemannian Lbb                & Constrained Optimization     \n");
fprintf("IT     & RES      & T(s)     & IT     & RES      & T(s)     & IT     & RES      & T(s)     \n");
i = 1;
%A = testmatrices{i};
load('../Matrices/GD96_c.mat')
A = Problem.A;
G = graph(A,'omitselfloops');
n = size(A,1);
D = spdiags(G.degree,0,n,n);
A = full(D\A);
[pi,~] = eigs(A.',1,'largestabs');
pi = pi./norm(pi,1);
pi = pi*sign(pi(1));
n = size(A,1);
optionsolve = initoptions();
optionsolve.method = "cg";
optionsolve.verbose = true;
optionsolve.correction = true;
optionsolve.formulation = "block";
optionsolve.plot = true;
optionsolve.storeiter = false;
M = multinomialfixedstochasticfactory(pi,optionsolve);
X0 = M.rand(); % Initial guess
% Riemannian optimization algorithm
problem.M = M;
problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p)-A);
problem = manoptAD(problem);
options.tolgradnorm = 1e-3;%tolerances(i);
options.verbosity = 2;
options.maxiter = 50;
[X1, xcost, info, options] = trustregions(problem,X0,options);
residual(1,i) = norm(mpower(X1,p)-A,"fro");
time(1,i) = info.time;
iter(1,i) = length([info.iter]);
