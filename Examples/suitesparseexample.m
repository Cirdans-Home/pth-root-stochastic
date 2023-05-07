%% Checking the stationary distribution when optimizing with MANOPT
% This example uses the Markov chain induced by a random walk on an
% undirected graph. We compare the stationary vector of the original matrix
% and the one we obtain for the Riemannian square root approximation on the
% multinomial stochastic manifold.
%
% To produce the figures for the manuscript we use the matlab2tikz
% function. It is put in a try/cat statement, so if you have not installed
% it, you can run everything as is and not have the figures.

clear; clc;

try manopt_version;
    fprintf("Using MANOPT version %d.%d.%d\n\n",manopt_version);
catch
    fprintf("Loading MANOPT\n");
    here = pwd;
    cd ../manopt
    addpath('auxiliaries')
    importmanopt;
    cd(here);
end

load('../Matrices/GD96_c.mat')
A = Problem.A;
G = graph(A,'omitselfloops');
n = size(A,1);
D = spdiags(G.degree,0,n,n);
A = D\A;
e = ones(n,1);

[pi,~] = eigs(A.',1,'largestabs');
pi = pi./norm(pi,1);

%% Building variety
manifold = multinomialfactory(n,n); % These are column stochastic!
%% Building the problem
p = 2;
A = full(A);
problem.M = manifold;
problem.cost = @(x) 0.5*cnormsqfro(mpower(x,p).'-A);
problem = manoptAD(problem);
options.tolgradnorm = 1e-7;
%% Solving
[x, xcost, info, ~] = trustregions(problem,[],options);


%% Convergence history
figure(1)
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

%% Perron vector
[pix,ev] = eigs(x,1,'largestabs');
pix = sign(pix(1))*pix./norm(pix,1);

figure(2)
plot(1:n,pix,'X',1:n,pi,'o')

%% Matrix A
figure(3)
subplot(3,2,1)
plot(G,'XData',Problem.aux.coord(:,1),...
    'YData',Problem.aux.coord(:,2))
axis square
set(gca,"Box","off","XTick",[],"YTick",[])
subplot(3,2,2)
spy(A)
subplot(3,2,3)
plot(G,'XData',Problem.aux.coord(:,1),...
    'YData',Problem.aux.coord(:,2),'MarkerSize',pi*5*n)
axis square
set(gca,"Box","off","XTick",[],"YTick",[])
title('Original Matrix')
subplot(3,2,4)
plot(G,'XData',Problem.aux.coord(:,1),...
    'YData',Problem.aux.coord(:,2),'MarkerSize',pix*5*n)
axis square
set(gca,"Box","off","XTick",[],"YTick",[])
title('Approximate Root')
subplot(3,2,[5,6])
plot(1:n,pi,'X',1:n,pix,'o')
xlim([1,n])
xlabel('node')
ylabel('Stationary Distribution')
legend('Original Matrix','Approximate root')
try
    % If export fig command is available
    set(gcf,'Color','none')
    export_fig('eigenfailure.pdf')
catch
    % Matlab's default
    print(gcf,'-dpdf','eigenfailure.pdf');
end

% Plot only the eigenvector
figure(4)
plot(1:n,pi,'X',1:n,pix,'o')
xlim([1,n])
xlabel('node')
ylabel('Stationary Distribution')
legend('Original Matrix','Approximate root')
try
    matlab2tikz('eigenvector.tikz')
catch
    warning("This output needs the matlab2tikz code.")
end