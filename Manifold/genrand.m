function A = genrand(pi)
%GENDRAND Generates a Random row-stochastic matrix with given left-hand
%side eigenvector pi
%   The function uses the decomposition of A as
%   A = e*pi' + B with B 1 = 0 and B < - e*pi' (componentwise) 
n = length(pi);
e = ones(n,1);
l = randn(n-1,1);
l = l./(2*norm(l,inf));
V = [e,null(e')].';
Q = [pi,null(pi')].';
t = 1;
A = e*pi' + t*Q'*diag([0;l])*V;
while any(A(:) < 0)
    t = t/2;
    A = e*pi' + t*Q'*diag([0;l])*V;
end
end

