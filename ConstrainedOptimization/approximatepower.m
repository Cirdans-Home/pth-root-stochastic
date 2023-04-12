function [X,output,history] = approximatepower(A,p,A0,tol,maxIterations)
%%APPROXIMATEPOWER solves the constrained optimization problem
%   X = arg min || X^p - A ||_F
% requiring X to be row-stochastic, i.e., X >= 0 and X 1 = 1.

n = size(A,1);
LB = zeros(n^2,1);
I = speye(n,n);
e = ones(n,1);
Aeq = kron(e',I); 
beq = ones(n,1);
history = StoreOptimizationValues(maxIterations+1,n^2);

options = optimset("Algorithm","interior-point",...
    "Display","none",...
    "TolFun",tol,...
    "OutputFcn",...
    @(x,optimValues,state) outfun(x,optimValues,state,history));
tic;
[X,~,~,output] = fmincon(@(X) objective(X,A,p),A0(:),...
    [],[],Aeq,beq,LB,[],[],options);
output.time = toc;
X = reshape(X,n,n);
% Remove unused space
history.targetFunctionValues = history.targetFunctionValues(1:output.iterations);
history.variableValues = history.variableValues(1:output.iterations,:);


end

function Y = objective(X,A,p)
%%OBJECTIVE is the objective function to minimize

n = sqrt(length(X));
Xp = mpower(reshape(X,n,n),p);
Y = 0.5*norm(Xp(:)-A(:),"fro")^2;

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