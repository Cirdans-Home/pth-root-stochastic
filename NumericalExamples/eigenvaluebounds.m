%% Check of the bound on the Eigenvalues
% Show the bound on the eigenvalues of the matrix for different
% manifolds/stationary vectors.

clear; clc; close all;

addpath("../Manifold/");
addpath("../Matrices/");
addpath("../Utilities/")

%% Generate test matrices
seed = 100;
classes = 1:6; % All types
numbers = ones(size(classes));
matrixset = matrixgenerator(100,2,seed,classes,numbers);

for i=1:length(classes)
    pi = pvgth(matrixset{i});

    % Build Manifold
    manifold = multinomialfixedstochasticfactory(pi);
    % Generate point on manifold
    A = manifold.rand();
    % Matrix of which we need to estimate the eigenvalues
    I = eye(size(A));
    Dpi = spdiags(pi,0,size(A,1),size(A,2));
    deltavec = A'*Dpi*pi;
    delta = spdiags( deltavec,0,size(A,1),size(A,2));
    M = [I, Dpi*A; A'*Dpi, delta];
    % Upper bound
    ubound = max([1+norm(pi,"inf"),2*norm(pi,"inf")]);
    % Lower bound
    [rstar,k] = min(max(1-A,[],1));
    deltastar = min(deltavec((1:size(A,1)~=k)));
    if rstar + deltastar < 1
        lbound = deltastar*(1 - rstar/(1-deltastar));
    else
        lbound = (deltastar+(1-sqrt(deltastar*(deltastar+4*rstar-2)+1)))/2;
    end
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

figure(1)
set(gca,'XScale','log')
set(gca,'XGrid','on')
yticks(1:length(classes))
yticklabels({'Unif.','pth power of unif.','exp of intensity','K80 (Emb.)','K80 (Not. Emb.)','Pei'})