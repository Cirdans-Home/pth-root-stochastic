function options = initoptions()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

options.formulation = "block"; % schur
options.method = "cg";         % lsqr svd 
options.correction = false;    % Do rank 1 update?
options.verbose = false;       % Prints info on linear system solution
options.plot = false;          % Plot convergence of lineas system solver
options.fighandle = [];        % New-handle o reuse?

end