%% Optimization of row-stochastic matrices
% Performance tests for the solution of the approximation problem on the
% set of row-stochastic matrices. Comparison between the usage of the
% Riemannian optimization algorithms and the standard constrained
% optimization one.

clear; clc; close all hidden;

try manopt_version;
    fprintf("Using MANOPT version %d.%d.%d\n\n",manopt_version);
catch
    fprintf("Loading MANOPT\n");
    here = pwd;
    cd ../manopt
    importmanopt;
    cd(here);
end
addpath('../Utilities/'); % code for performance profile plots
addpath('../Matrices/'); % code for the generation of test matrices
addpath('../ConstrainedOptimization/'); % code for constrained optimization
%% Generation of the test problems
n    = 100; % Size of the test matrices: some are of fixed size!
seed = 42;  % Random number generator seed use it to enforce reproducibility
p = 2;      % Square-root, change it to look for other roots
classes = 1:6;
numtestperclass = 40;
numbers = numtestperclass*ones(length(classes),1);
testmatrices = matrixgenerator(n,p,seed,classes,numbers);
tolerances = kron([1e-3,1e-6,1e-6,1e-6,1e-3,1e-6].',ones(numtestperclass,1));

%% Matrices for the collection of the performances
time = zeros(3,length(testmatrices));
iter = zeros(3,length(testmatrices));
residual = zeros(3,length(testmatrices));
%% Running optimization tests
fprintf("Riemannian trust region      & Riemannian Lbb                & Constrained Optimization     \n");
fprintf("IT     & RES      & T(s)     & IT     & RES      & T(s)     & IT     & RES      & T(s)     \n");
advbar = waitbar(0,"Running tests");
for i=1:length(testmatrices) % Test of all the matrices
    A = testmatrices{i};
    n = size(A,1);
    M = multinomialfactory(n,n);
    X0 = M.rand(); % Initial guess
    % Riemannian optimization algorithm
    problem.M = M;
    problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
    problem = manoptAD(problem);
    options.tolcost = 1e-3;
    %options.tolgradnorm = 1e-3;%tolerances(i);
    options.verbosity = 0;
    [X1, xcost, info, options] = trustregions(problem,X0,options);
    residual(1,i) = norm(mpower(X1.',p)-A,"fro");
    time(1,i) = info.time;
    iter(1,i) = length([info.iter]);
    clear problem options
    % Riemannian optimization algorithm
    problem.M = M;
    problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
    problem = manoptAD(problem);
    options.tolcost = 1e-3;
    %options.tolgradnorm = 1e-3;%tolerances(i);
    options.verbosity = 0;
    options.strategy = 'alternate';
    [X3, xcost3, info, options] = rlbfgs(problem,X0,options);
    residual(2,i) = norm(mpower(X1.',p)-A,"fro");
    time(2,i) = info.time;
    iter(2,i) = length([info.iter]);
    clear problem M
    % Constrained optimization algorithm
    [X2,output,history] = approximatepower(A,p,X0,1e-3,1000); % tolerances(i)
    residual(3,i) = norm(mpower(X2,p)-A,"fro");
    time(3,i) = output.time;
    iter(3,i) = output.funcCount;
    % Table
    fprintf('%06d & %1.2e & %1.2e & %06d & %1.2e & %1.2e & %06d & %1.2e & %1.2e \n',...
        iter(1,i),residual(1,i),time(1,i),...
        iter(2,i),residual(2,i),time(2,i),...
        iter(3,i),residual(3,i),time(3,i))
    waitbar(i/length(testmatrices),advbar,"Running tests");
end
close(advbar);

%% Plot Performance profiles
figure(1);
perfprof(time.');
legend({'Riemannian trust-region','Riemannian LBFGS', ...
    'Interior point method'},'Location','southeast')
title('Execution time (s)')
figure(2);
perfprof(residual.');
legend({'Riemannian trust-region','Riemannian LBFGS', ...
    'Interior point method'},'Location','southeast')
title('Residual')