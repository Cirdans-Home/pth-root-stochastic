function M = multinomialfixedstochasticfactory(pi,optionsolve)
% Manifold of n-by-n row stochastic matrices with positive entries and
% fixed left-hand stochastic eigenvector
%
% function M = multinomialfixedstochasticfactory(n,pi)
%
% M is a Manopt manifold structure to optimize over the set of n-by-n
% matrices with (strictly) positive entries and such that the entries of
% each row sum to one and the left eigenvector is pi.
%
% Points on the manifold and tangent vectors are represented naturally as
% matrices of size n. The Riemannian metric imposed on the manifold is the
% Fisher metric, that is, if X is a point on the manifold and U, V are two
% tangent vectors:
%
%     M.inner(X, U, V) = <U, V>_X = sum(sum(U.*V./X)).
%
% The retraction here provided is only first order. Consequently, the
% slope test in the checkhessian tool is only valid at points X where the
% gradient is zero. Furthermore, if some entries of X are very close to
% zero, this may cause numerical difficulties that can also lead to a
% failed slope test. More generally, it is important that the solution of
% the optimization problem should have positive entries, sufficiently far
% away from zero to avoid numerical issues.
%
% This file is based on the multinomialdoublystochasticfactory.m code by
% Ahmed Douik and Nicolas Boumal.
%

    n = length(pi);
    e = ones(n, 1);

    maxDSiters = 100 + 2*n;
    
    function [alpha, beta] = mylinearsolve(X, b)
        % zeta = sparse(A)\b; % sparse might not be better perf.-wise.
        % where A = [eye(n) X ; X' eye(n)];
        %
        % For large n use the pcg solver instead of \.
        % [zeta, ~, ~, iter] = pcg(sparse(A), b, 1e-6, 2*n);
        %
        % Even faster is to create a function handle
        % computing A*x (x is a given vector). 
        % Make sure that A is not created, and X is only 
        % passed with mylinearsolve and not A.

        switch upper(optionsolve.formulation)
            case "BLOCK"
                % We solve here the system in its 2 x 2 block formulation
                switch upper(optionsolve.method)
                    case "CG"
                        [zeta, flag, res, iter,resvec] = ...
                            pcg(@(x) mycompute(x,true), b, 1e-6, 2*n);
                    case "LSQR"
                        [zeta, flag, res, iter,resvec] = ...
                            lsqr(@mycompute, b, 1e-6, 2*n);
                    case "SVD"
                        S1 = [eye(size(X)) diag(pi)*X;
                            X.'*diag(pi) diag(X.'*(diag(pi)*pi))];
                        zeta = pinv(S1)*b;
                    otherwise
                        error("Manopt:unknown_linear_solver")
                end
            case "SCHUR"
                bschur = b(n+1:end,1) - X.'*diag(pi)*b(1:n,1);
                % We solve here the system in its reduced formulation
                switch upper(optionsolve.method)
                    case "CG"
                        [zetaschur, flag, res, iter,resvec] = ...
                            pcg(@(x) mycompute_schur(x,true), bschur, ...
                            1e-6, 2*n); 
                    case "LSQR"
                        [zetaschur, flag, res, iter,resvec] = ...
                            lsqr(@mycompute_schur, bschur, ...
                            1e-6, 2*n);
                    case "SVD"
                        S1schur = diag(X.'*(diag(pi)*pi)) - X.'*diag(pi.^2)*X;
                        zetaschur = pinv(S1schur)*bschur;
                    otherwise
                        error("Manopt:unknown_linear_solver")
                end
                % reconstruct the whole vector
                zeta = [b(1:n,1) - pi.*(X*zetaschur);...
                    zetaschur];
            otherwise
                error("Manopt:unknown_linear_formulation")
        end

        if (flag > 0 || optionsolve.verbose)
            fprintf("%s & %s & %d & %1.2e \\\\\n",optionsolve.formulation,...
                optionsolve.method,iter,res);
            if (optionsolve.plot)
                if isempty(optionsolve.fighandle)
                    optionsolve.fighandle = figure();
                else 
                    figure(optionsolve.fighandle);
                    hold on
                end
                if optionsolve.correction
                    correctionstring = '+ rank 1';
                else 
                    correctionstring = '';
                end
                semilogy(1:length(resvec),resvec,'LineWidth',2,'DisplayName',...
                    sprintf('%s-%s %s',optionsolve.method, ...
                    optionsolve.formulation,correctionstring))
                hold off
            end
        end


        function Ax = mycompute(x,flag)
            % The matrix is symmetric, thus the flag is an empty
            % placeholder for the LSQR request for the transposition
            xtop = x(1:n,1);
            xbottom = x(n+1:end,1);
            Axtop = xtop + diag(pi)*X*xbottom;
            Axbottom = X'*diag(pi)*xtop + diag(pi'*diag(pi)*X)*xbottom;
            Ax = [Axtop; Axbottom];
            if optionsolve.correction
                Ax = Ax + ((pi'*xtop - e'*xbottom)/(pi'*pi + n))*[pi;-e];
            end
        end

        function Ax = mycompute_schur(x,flag)
            % The matrix is symmetric, thus the flag is an empty
            % placeholder for the LSQR request for the transposition
            Ax = (X.'*(diag(pi)*pi)).*x - X.'*(pi.^2.*(X*x));  
            if optionsolve.correction
                Ax = Ax + (sum(x)/(n))*e;
            end
        end
        
        alpha = zeta(1:n, 1);
        beta = zeta(n+1:end, 1);
    end

    M.name = @() sprintf(['%dx%d row-stochastic matrices with positive ' ...
        'entries and fixed left-stochastic eigenvector'], n, n);

    M.dim = @() (n-1)^2;

    % Fisher metric
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        ip = sum((eta(:).*zeta(:))./X(:));
    end

    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));

    M.dist = @(X, Y) error(['multinomialfixedstochasticfactory.dist not ' ...
        'implemented yet.']);

    % The manifold is not compact as a result of the choice of the metric,
    % thus any choice here is arbitrary. This is notably used to pick
    % default values of initial and maximal trust-region radius in the
    % trustregions solver.
    M.typicaldist = @() n;

    % Pick a random point on the manifold
    M.rand = @random;
    function X = random()
        X = abs(randn(n, n));
        X = modifiedsinkhorn(X,pi,maxDSiters);
    end

    % Pick a random vector in the tangent space at X.
    M.randvec = @randomvec;
    function eta = randomvec(X) % A random vector in the tangent space
        % A random vector in the ambient space
        Z = randn(n, n);
        % Projection of the vector onto the tangent space
        b = [sum(Z, 2) ; Z'*pi];
        [alpha, beta] = mylinearsolve(X, b);
        eta = Z - (alpha*e' + pi*beta').*X;
        % Normalizing the vector
        nrm = M.norm(X, eta);
        eta = eta / nrm;
    end

    % Projection of vector eta in the ambient space to the tangent space.
    M.proj = @projection; 
    function etaproj = projection(X, eta) % Projection of the vector eta in the ambeint space onto the tangent space
        b = [sum(eta, 2) ; eta'*pi];
        [alpha, beta] = mylinearsolve(X, b);
        etaproj = eta - (alpha*e' + pi*beta').*X;
    end

    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;

    % Conversion of Euclidean to Riemannian gradient
    M.egrad2rgrad = @egrad2rgrad;
    function rgrad = egrad2rgrad(X, egrad) % projection of the euclidean gradient
        mu = (X.*egrad); 
        b = [sum(mu, 2) ; mu'*pi];
        [alpha, beta] = mylinearsolve(X, b);
        rgrad = mu - (alpha*e' + pi*beta').*X;
    end

    % First-order retraction
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y = X.*exp(t*(eta./X));
        Y = modifiedsinkhorn(Y,pi);
        Y = max(Y, eps);
    end

    % Conversion of Euclidean to Riemannian Hessian
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, eta)
        
        % Computing the directional derivative of the Riemannian
        % gradient
        gamma = egrad.*X;
        gammadot = ehess.*X + egrad.*eta;
        
        bdot = [ gammadot*e ; gammadot.'*pi];
        b = [gamma*e ; gamma.'*pi];
        
        [alpha, beta] = mylinearsolve(X, b);
        S1 = [ zeros(size(eta)) , diag(pi)*eta; eta'*diag(pi) , diag(eta'*diag(pi)*pi) ];
        %S2 = [eye(size(X)) , diag(pi)*X; X'*diag(pi), diag( X'*diag(pi)*pi )];
        %[C,R] = qr(sparse([zeros(size(S1)), S2; S2, S1]),[bdot;b]);
        %LL = chol(S2);
        %sol = pinv([zeros(size(S1)), S2; S2, S1])*[b;bdot];
        %sol = lsqr([zeros(size(S1)), S2; S2, S1],[b;bdot],1e-6,100);
        %alphadot = sol(1:n);
        %betadot = sol(n+1:2*n);
        %alpha = sol(2*n+1:3*n);
        %beta = sol(3*n+1:end);
       
        
        [alphadot, betadot] = mylinearsolve(X, bdot - S1*[alpha; beta]); % %- [eta*beta; eta'*alpha]
        
        S = (alpha*e' + pi*beta');
        deltadot = gammadot - (alphadot*e' + pi*betadot').*X- S.*eta; % rgraddot

        % Computing Riemannian gradient
        delta = gamma - S.*X; % rgrad

        % Riemannian Hessian in the ambient space
        nabla = deltadot - 0.5*(delta.*eta)./X;

        % Riemannian Hessian on the tangent space
        rhess = projection(X, nabla);
    end


    % Miscellaneous manifold functions
    M.hash = @(X) ['z' hashmd5(X(:))];
    M.lincomb = @matrixlincomb;
    M.zerovec = @(X) zeros(n, n);
    M.transp = @(X1, X2, d) projection(X2, d);
    M.vec = @(X, U) U(:);
    M.mat = @(X, u) reshape(u, n, n);
    M.vecmatareisometries = @() false;
    
end
