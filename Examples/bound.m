%% NUMERICAL EXAMPLES OF THE BOUND FOR THE PROJECTION LINEAR SYSTEM
% Numerical examples for the bound in Proposition 2.

clc; clear; close all;



%load('../Matrices/gre_115.mat')
load('../Matrices/GD96_c.mat')
A = spones(Problem.A);
G = digraph(A,'omitselfloops');
n = size(A,1);
e = ones(n,1);
D = spdiags(A*e,0,n,n);
A = D\A;
I = speye(size(A));


[pi,~] = eigs(A.',1,'largestabs');
pi = pi./norm(pi,1);

addpath('../Manifold/')
manifold = multinomialfixedstochasticfactory(pi);

for i=1:20
    A = manifold.rand();

    Dpi = spdiags(pi,0,size(A,1),size(A,2));
    deltavec = A'*Dpi*pi;
    delta = spdiags( deltavec,0,size(A,1),size(A,2));
    M = [I, Dpi*A; A'*Dpi, delta];

    % Upper bound
    ubound = max([1+norm(pi,"inf"),2*norm(pi,"inf")]);
    % Lower bound
    [rstar,k] = max(min(A,[],1));
    deltastar = min(deltavec((1:size(A,1)~=k)));
    %if rstar + deltastar <= 1
    %    lbound = min(pi);
    %else
    lbound = (deltastar+(1-sqrt(deltastar*(deltastar+4*rstar-2)+1)))/2;
    %end
    % Eigenvalues
    evM = eig(full(M),'balance');
    evM(1) = 0;
    yone = ones(size(M,1),1);

    figure(1)
    hold on
    semilogx(evM,i*yone,'k.',ubound,i,'r|',lbound,i,'b|','MarkerSize',10)
    hold off

    fprintf('%e & %e & %d & %e & %e & %d \n', ...
        min(evM(evM > 0)),lbound,min(evM(evM > 0)) > lbound, ...
        max(evM),ubound,max(evM)<ubound);
    difflower(i) = (min(evM(evM > 0)) - lbound)/min(evM(evM > 0));
    diffupper(i) = (ubound - max(evM))/max(evM);
end

%% Plot of the bounds
figure(1)
set(gca,"XScale","log")
grid on

%% It the bound tight?
figure(2)
semilogy(1:length(difflower),difflower,'kx-',...
    1:length(diffupper),diffupper,'r.-','LineWidth',2)
legend('Lower bound','Upper bound')
title('Sharpness of the bound')
axis tight
ylim([1e-4,1e-1])

%% Rank 1 cases

A = e*pi';
Dpi = spdiags(pi,0,size(A,1),size(A,2));
deltavec = A'*Dpi*pi;
delta = spdiags( deltavec,0,size(A,1),size(A,2));
M = [I, Dpi*A; A'*Dpi, delta];

% Upper bound
ubound = max([1+norm(pi,"inf"),2*norm(pi,"inf")]);
% Lower bound
[rstar,k] = max(min(A,[],1));
deltastar = min(deltavec((1:size(A,1)~=k)));
%if rstar + deltastar <= 1
%    lbound = min(pi);
%else
lbound = (deltastar+(1-sqrt(deltastar*(deltastar+4*rstar-2)+1)))/2;
%end
% Eigenvalues
evM = eig(full(M),'balance');
evM(1) = 0;
yone = ones(size(M,1),1);

figure(3)
i = 1;
semilogx(evM,i*yone,'k.',ubound,i,'r|',lbound,i,'b|','MarkerSize',10)
ylim([0.8 1.2])
axis tight
