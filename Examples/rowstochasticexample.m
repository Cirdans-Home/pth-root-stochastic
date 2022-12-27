%% MARKOV INTERMEDIATE STEP FROM RIEMANNIAN GEOMETRY

clear; clc;

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
a = 1/5;
A = 1/3*[1-2*a 1+a 1+a
   1+a 1-2*a 1+a
   1+a 1+a 1-2*a];
n = size(A,1);
p = 2;
e = ones(n,1);

%% Building variety
manifold = multinomialfactory(n,n); % These are column stochastic!
%% Building the problem
problem.M = manifold;
problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
problem = manoptAD(problem);
options.tolgradnorm = 5e-9;
[x, xcost, info, options] = trustregions(problem,[],options);

%% Convergence history
figure(1)
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');