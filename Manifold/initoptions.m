function options = initoptions()
%INITOPTIONS Setup options for the linear solvers inside the multinomial
%fixed stochastic factory manifold. It permits to select in between the
%different formulation of the linear system needed for the projection, and
%to store the relevant quantities for convergence analysis purposes.

options.formulation = "block"; % schur
options.method = "direct";     % cg pcg1 pcg2 lsqr direct 
options.threshold = 1e-2;      % Threshold for both MIC and dropping in the 
                               % Neumann series preconditioner
options.kappa = 3;             % Number of terms in the Neumann series
options.verbose = false;       % Prints info on linear system solution
options.plot = false;          % Plot convergence of lineas system solver
options.fighandle = [];        % New-handle o reuse?
options.storeiter = false;     % Do I have to store linear iteration history?
                               % This activates a *very slow* procedure!
                               % Use only for debug/study purposes

end