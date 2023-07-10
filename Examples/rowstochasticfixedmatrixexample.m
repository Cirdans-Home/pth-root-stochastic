%% MARKOV INTERMEDIATE STEP FROM RIEMANNIAN GEOMETRY

clear; clc; close all;

try manopt_version;
    fprintf("Using MANOPT version %d.%d.%d\n\n",manopt_version);
catch
    fprintf("Loading MANOPT\n");
    here = pwd;
    cd ../manopt
    importmanopt;
    cd(here);
end

%% We start from the Markov function we are interested in
a = 1/6;
A = 1/3*[1-2*a 1+a 1+a
   1+a 1-2*a 1+a
   (1+a)/3 (1+a) 3-(1+a)-(1+a)/3];
n = size(A,1);
p = 2;
e = ones(n,1);

A0 = rand(size(A));
D = diag(sum(A0,2));
A0 = D\A0;

%% Building manifold
manifold = multinomialfactory(n,n); % These are column stochastic!
%% Building the problem
problem.M = manifold;
problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
problem = manoptAD(problem);
options.tolgradnorm = 5e-9;
[X1, xcost, info, options] = trustregions(problem,A0,options);
fprintf("(Riemannian) ||X 1 - 1 || = %e\n",norm(X1.'*e-e,"inf"));
fprintf("||X^p - A|| = %e\n",norm(mpower(X1,p).'-A,"fro"));

%% Constrained optimization solution
[X2,output,history] = approximatepower(A,p,A0,max(xcost,eps),1000);
fprintf("(fmincon) ||X 1 - 1 || = %e\n",norm(X2*e-e,"inf"));
fprintf("||X^p - A|| = %e\n",norm(mpower(X2,p)-A,"fro"));

%% Constrained optimization solution with fixed stationary state
clear options problem manifold
pi = pvgth(A);
optionsolve = initoptions();
optionsolve.formulation = "schur";
optionsolve.method = "cg";
optionsolve.correction = true;
optionsolve.verbose = true;
manifold2 = multinomialfixedstochasticfactory(pi,optionsolve);
problem2.M = manifold2;
problem2.cost = @(x) 0.5*cnormsqfro(mpower(x,p)-A);
problem2 = manoptAD(problem2);
options.tolgradnorm = 5e-9;
[X3, xcost, info, options] = trustregions(problem2,A0,options);
fprintf("(Riemannian) ||X 1 - 1 || = %e\n",norm(X3*e-e,"inf"));
fprintf("||X^p - A|| = %e\n",norm(mpower(X3,p).'-A,"fro"));

%% Convergence history
figure(1)
semilogy([info.iter]+1, [info.cost], '.-',...
    1:output.iterations,history.targetFunctionValues,'--',...
    'LineWidth',2);
xlabel('Iteration number');
ylabel('Objective function value')
legend('Manifold Trust-Region','Interior-Point (Constrained Optimization)');
axis tight
axis square
grid on