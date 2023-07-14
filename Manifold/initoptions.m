function options = initoptions()
%INITOPTIONS Setup options for the linear solvers inside the multinomial
%fixed stochastic factory manifold. It permits to select in between the
%different formulation of the linear system needed for the projection, and
%to store the relevant quantities for convergence analysis purposes.

options.formulation = "block"; % schur
options.method = "svd";        % cg lsqr svd 
options.correction = false;    % Do rank 1 update?
options.verbose = false;       % Prints info on linear system solution
options.plot = false;          % Plot convergence of lineas system solver
options.fighandle = [];        % New-handle o reuse?
options.storeiter = false;     % Do I have to store linear iteration history?
                               % This activates a *very slow* procedure!
                               % Use only for debug/study purposes

end