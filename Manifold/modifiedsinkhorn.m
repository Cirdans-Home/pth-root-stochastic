function [B,u,v] = modifiedsinkhorn(A,pi,maxit,checkperiod)
%%MODIFIEDSINKHORN applies the modified Sinkhorn algorithm to get a
%%row-stochastic matrix with stationary distribution pi.

N = size(A,1);
tol = eps(N);

% Set default values
if ~exist('maxit', 'var') || isempty(maxit)
    maxit = 1000;
end
if ~exist('checkperiod', 'var') || isempty(checkperiod)
    checkperiod = round(N/10);
end

Ahat = diag(pi)*A;
% Number of iteration
iter = 0;
% Initialize u and v
u = ones(N,1); v = ones(N,1); e = v;
while iter < maxit
    iter = iter + 1;

    % previous for safekeeping
    u_prev = u;
    v_prev = v;

    row = Ahat*v;

    % update u and v
    u = pi./(row);
    v = pi./(Ahat'*u);

    % Check if converged
    if mod(iter, checkperiod) == 0
        gap = abs(u'*row - 1);
        if isnan(gap)
            break;
        end
        if gap <= tol
            if norm(diag(u)*A*diag(v)*e-e,"inf") < tol && ...
                    norm(pi'*diag(u)*A*diag(v)-pi',"inf") < tol
                break;
            end
        end
    end

    if any(isinf(u)) || any(isnan(u)) || any(isinf(v)) || any(isnan(v))
        warning('DoublyStochasticProjection:NanInfEncountered', ...
            'Nan or Inf occured at iter %d. \n', iter);
        u = u_prev;
        v = v_prev;
        break;
    end


end

% The matrix we want is built from A as in the Theorem
B = diag(u)*A*diag(v);


end