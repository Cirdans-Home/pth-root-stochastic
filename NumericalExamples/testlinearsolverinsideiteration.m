%% Optimization with stochastic matrices with assigned pi vector
% Performance tests for the solution of the approximation problem on the
% set of row-stochastic matrices with assigned stationary vector. Comparison
% between the usage of the Riemannian optimization algorithms and the
% standard constrained optimization one.

clear all; clc;

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
rng(seed);  % Fix the randoom number generator
%% Select a matrix
load('../Matrices/GD96_c.mat')
%load('../Matrices/gre_115.mat')
%load('../Matrices/Roget.mat')
%% Build the random walk on the graph
A = Problem.A;
G = digraph(A);
[bin,binsize] = conncomp(G,'Type','strong');
B = binsize(bin) == max(binsize);
n = size(A,1);
A = A(B,B);
fprintf('%s & %s & %d & ',Problem.name,Problem.kind,n);
n = size(A,1);
D = spdiags(A*ones(n,1),0,n,n);
A = full(D\A);
ev = eig(full(A));
[pi,~] = eigs(A.',1,'largestabs');
pi = pi./norm(pi,1);
pi = pi*sign(pi(1));
fprintf('%d \\\\\n',n);
optionsolve = initoptions();
optionsolve.method = "pcg2";
optionsolve.verbose = true;
optionsolve.threshold = 1e-2;
optionsolve.kappa = 3;
optionsolve.formulation = "schur";
optionsolve.fighandle = [];
optionsolve.plot = false;
optionsolve.storeiter = true;
M = multinomialfixedstochasticfactory(pi,optionsolve);
X0 = M.rand(); % Initial guess
%% Riemannian optimization algorithm
problem.M = M;
problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p)-A);
problem = manoptAD(problem);
options.tolgradnorm = 1e-3;
options.verbosity = 2;
options.maxiter = 50;
[X1, xcost, info, options] = trustregions(problem,X0,options);
residual = norm(mpower(X1,p)-A,"fro");
time = info.time;
iter = length([info.iter]);
