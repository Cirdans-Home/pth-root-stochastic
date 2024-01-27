function [X,output,history] = approximatepower(A,p,A0,tol,maxIterations,varargin)
%%APPROXIMATEPOWER solves the constrained optimization problem
%   X = arg min || X^p - A ||_F
% requiring X to be row-stochastic, i.e., X >= 0 and X 1 = 1.
% The optional argument can be either 'FD' to use finite difference
% approximations of Hessian and Gradients, or 'ANALYTICAL' to use the one
% computed.

n = size(A,1);
LB = zeros(n^2,1);
I = speye(n,n);
e = ones(n,1);
Aeq = kron(e',I);
beq = ones(n,1);
history = StoreOptimizationValues(maxIterations+1,n^2);

if nargin > 5
    whathessian = varargin{1};
else
    whathessian = "FD";
end

switch upper(whathessian)
    case "FD"
        options = optimset('Algorithm','interior-point',...
            'Display','none',...
            "TolFun",tol,...
            "OutputFcn",...
            @(x,optimValues,state) outfun(x,optimValues,state,history));
    case "ANALYTICAL"
        HessMultFcn = @(x,lambda,v) hessian(x,lambda,v,A,p);
        options = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','none', ...
            'SpecifyObjectiveGradient',true, ...
            'SpecifyConstraintGradient',false, ...
            'SubproblemAlgorithm','cg', ...
            'TolFun',tol,...
            'HessianMultiplyFcn', @(x,lambda,v) HessMultFcn(x,lambda,v), ...
            'OutputFcn',...
            @(x,optimValues,state) outfun(x,optimValues,state,history));
    otherwise
        warning("Defaulting to Finite Differences")
        options = optimset('fmincon',...
            "Algorithm","interior-point",...
            "Display","none",...
            'SpecifyObjectiveGradient',true,...
            'SpecifyConstraintGradient',true,...
            'SubproblemAlgorithm','cg',...
            'HessianMultiplyFcn',@HessMultFcn,...
            "TolFun",tol,...
            "OutputFcn",...
            @(x,optimValues,state) outfun(x,optimValues,state,history));
end



tic;
[X,~,~,output] = fmincon(@(X) objective(X,A,p),A0(:),...
    [],[],Aeq,beq,LB,[],[],options);
output.time = toc;
X = reshape(X,n,n);
% Remove unused space
history.targetFunctionValues = history.targetFunctionValues(1:output.iterations);
history.variableValues = history.variableValues(1:output.iterations,:);


end

function [Y,g] = objective(X,A,p)
%%OBJECTIVE is the objective function to minimize

n = sqrt(length(X));
Xm = reshape(X,n,n);
Xp = mpower(Xm,p);
Y = 0.5*norm(Xp(:)-A(:),"fro")^2;

if nargout > 1 % Gradient is required
    F = Xp - A;
    g = zeros(n,n);
    for j=1:p
        g = g + mpower(Xm',j-1)*F*mpower(Xm',p-j);
    end
end

end

function stop = outfun(x, optimValues, state, store)
stop = false;
switch state
    case 'iter'
        %Store important values in array
        i=optimValues.iteration+1;
        store.addValues(i,optimValues.fval,x);
end
end

function y = hessian(x,lambda,v,A,p)
%HESSIAN implementation, this uses the matrix formulation and a number of
%conversion between matrix and vector arguments

n = size(A,1);      % Size of the problem

X = reshape(x,n,n); % Current point
E = reshape(v,n,n); % Evaluation matrix

% Hessian relative to free objective
F = mpower(X,p)-A;
Y = zeros(n,n);
for j=1:p
    S = zeros(n,n);
    for l=1:p-j
        S = S+mpower(X',p-j-l)*E'*mpower(X',l-1);
    end
    Y = Y + mpower(X',j-1)*F*S;
    % Second term
    S = zeros(n,n);
    for k=1:p
        S = S + mpower(X,p-k)*E*mpower(X,k-1)*mpower(X',p-j);
    end
    Y = Y + mpower(X',j-1)*S;
    for i=1:j-1
        Y = Y + mpower(X',j-1-i)*E'*mpower(X',i-1)*F*mpower(X',p-j);
    end
end
% Here goes the part with the constraints that is always empty in our case.
l1 = lambda.ineqnonlin;
l2 = lambda.eqnonlin;

% Final Conversion
y = Y(:);

end