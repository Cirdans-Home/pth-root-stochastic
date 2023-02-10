function [B,u,v] = modifiedsinkhorn(A,pi,maxit,checkperiod)
%%MODIFIEDSINKHORN applies the modified Sinkhorn algorithm to get a
%%row-stochastic matrix with stationary distribution pi.

N = size(A,1);
tol = eps(N);

% Set default values
if ~exist('maxiter', 'var') || isempty(maxit)
    maxit = N^2;
end
if ~exist('checkperiod', 'var') || isempty(checkperiod)
    checkperiod = 100;
end

Ahat = diag(pi)*A;
% Number of iteration
iter = 0;
% Initialize u and v
u = ones(N,1); v = ones(N,1);
while iter < maxit
    iter = iter + 1;
    
    row = Ahat*v;
    % Check gap condition only at checkperiod intervals.
    % It saves computations for large-scale scenarios.
    if mod(iter, checkperiod) == 0
        gap = max(abs(row .* u - 1));
        if isnan(gap)
            break;
        end
        if gap <= tol
            break;
        end
    end
    % store old u and v in case of nasty breakdowns
    uprev = u;
    vprev = v;
    % update u and v
    u = pi./row;
    v = pi./(Ahat'*u);

    if any(isinf(v)) || any(isnan(v)) || any(isinf(u)) || any(isnan(u))
        warning('FixedStochasticProjection:NanInfEncountered', ...
            'Nan or Inf occured at iter %d. \n', iter);
        u = uprev;
        v = vprev;
        break;
    end

end

% The matrix we want is built from A as in the Theorem
B = diag(v)*A*diag(u);

end