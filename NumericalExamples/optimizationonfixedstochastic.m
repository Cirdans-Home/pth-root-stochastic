%% Optimization with stochastic matrices with assigned pi vector
% Performance tests for the solution of the approximation problem on the
% set of row-stochastic matrices with assigned stationary vector. Comparison
% between the usage of the Riemannian optimization algorithms and the
% standard constrained optimization one.

clear; clc; close all hidden;

try manopt_version;
    fprintf("Using MANOPT version %d.%d.%d\n\n",manopt_version);
catch
    fprintf("Loading MANOPT\n");
    here = pwd;
    cd ../manopt
    importmanopt;
    addpath("auxiliaries");
    cd(here);
end
addpath('../Utilities/'); % code for performance profile plots
addpath('../Matrices/'); % code for the generation of test matrices
addpath('../ConstrainedOptimization/'); % code for constrained optimization
addpath('../Manifold'); % code for the new manifold
%% Generation of the test problems
n    = 30; % Size of the test matrices: some are of fixed size!
seed = 42;  % Random number generator seed use it to enforce reproducibility
p = 2;      % Square-root, change it to look for other roots
classes = 1:6;
numtestperclass = 1;
numbers = numtestperclass*ones(length(classes),1);
testmatrices = matrixgenerator(n,p,seed,classes,numbers);
nclasses = length(classes);

%% Matrices for the collection of the performances
time = zeros(2,length(testmatrices));
iter = zeros(2,length(testmatrices));
residual = zeros(2,length(testmatrices));
residualpi = zeros(2,length(testmatrices));
%% Running optimization tests
fprintf("Riemannian trust region      & Riemannian Lbb                & Constrained Optimization     \n");
fprintf("IT     & RES      & T(s)     & IT     & RES      & T(s)     & IT     & RES      & T(s)     \n");
advbar = waitbar(0,"Running tests");
for i=1:length(testmatrices) % Test of all the matrices
    A = testmatrices{i};
    [pi,~] = eigs(A.',1,'largestabs');
    pi = pi./norm(pi,1);
    pi = pi*sign(pi(1));
    n = size(A,1);
    optionsolve = initoptions();
    optionsolve.method = "svd";
    optionsolve.verbose = false;
    M = multinomialfixedstochasticfactory(pi,optionsolve);
    M2 = multinomialfactory(n,n);
    X0 = M.rand(); % Initial guess
    % Riemannian optimization algorithm
    problem.M = M;
    problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p)-A);
    problem = manoptAD(problem);
    options.tolgradnorm = 1e-4;
    options.verbosity = 0;
    [X1, xcost, info, options] = trustregions(problem,X0,options);
    residual(1,i) = norm(mpower(X1,p)-A,"fro");
    time(1,i) = info.time;
    iter(1,i) = length([info.iter]);
    evX1 = pvgth(X1);
    residualpi(1,i) = norm(evX1-pi,"inf");
    clear problem options
    % Riemannian optimization algorithm
    problem.M = M2;
    problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
    problem = manoptAD(problem);
    options.tolgradnorm = 1e-4;
    options.verbosity = 0;
    options.strategy = 'alternate';
    [X2, xcost2, info, options] = trustregions(problem,X0,options);
    residual(2,i) = norm(mpower(X2.',p)-A,"fro");
    time(2,i) = info.time;
    iter(2,i) = length([info.iter]);
    evX2 = pvgth(X2);
    residualpi(2,i) = norm(evX2-pi,"inf");
    clear problem M
    % Table
    fprintf('%06d & %1.2e & %1.2e & %06d & %1.2e & %1.2e \n',...
        iter(1,i),residual(1,i),time(1,i),...
        iter(2,i),residual(2,i),time(2,i))
    waitbar(i/length(testmatrices),advbar,"Running tests");
end
close(advbar);

%% Residual on the stationary vector
hfig = figure('Position',[2839 844 1132 396]);
subplot(1,2,1)
semilogy(1:nclasses,residual(1,:),'x',1:nclasses,residual(2,:),'o','MarkerSize',7)
ylabel('$\| X^p - A \|_F$','Interpreter','latex')
legend({'$S_n^\pi$','$S_n$'},'Interpreter','latex','Location','eastoutside')
xticks(1:length(classes))
xticklabels({'Unif.','pth power of unif.','exp of intensity','K80 (Emb.)','K80 (Not. Emb.)','Pei'})
xtickangle(45)
axis tight
subplot(1,2,2)
semilogy(1:nclasses,residualpi(1,:),'x',1:nclasses,residualpi(2,:),'o','MarkerSize',7)
%xlabel('Matrix class')
ylabel({'Infinity norm error between';'target steady state and obtained one'},'Interpreter','latex')
legend({'$S_n^\pi$','$S_n$'},'Interpreter','latex','Location','eastoutside')
hold on
semilogy(linspace(1,nclasses,100),eps*ones(100,1),'k--','DisplayName','$\epsilon$')
hold off
axis tight
xticks(1:length(classes))
xticklabels({'Unif.','pth power of unif.','exp of intensity','K80 (Emb.)','K80 (Not. Emb.)','Pei'})
xtickangle(45)
%% Save figures to file
% If matlab2tikz is available produces tikz figures, otherwise it exports
% them as eps, see:
% https://github.com/matlab2tikz/matlab2tikz
try
    figure(hfig)
    matlab2tikz('filename','preserving-pi.tikz',...
        'height','1.5in',...
        'width','0.8\columnwidth', ...
        'parseStrings',false);
catch
    figure(1);
    savefig(gcf,'preserving-pi','epsc')
end
