function pv=pvgth(P)
%%PVGTH Computes the steady state vector x of the stochastic matrix P
% through LU factorization of I-P, with diagonal adjustment
n=size(P,1);
A=eye(n)-P;
% diagonal adjustement
for i=1:n
     A(i,i)=-sum(A(i,1:i-1))-sum(A(i,i+1:n));
end
% Gaussian elimination
for k=1:n-1
% Update L
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
% Update U
    for i=k+1:n
        A(i,k+1:n)=A(i,k+1:n)-A(i,k)*A(k,k+1:n);
    end
% Diagonal adjustment
    for i=k+1:n
        A(i,i)=-sum(A(i,k+1:i-1))-sum(A(i,i+1:n));
    end
end
% Back substitution
pv=ones(n,1);
for j=n-1:-1:1
    s=0;
    for i=j+1:n
        s=s+pv(i)*A(i,j);
    end
    pv(j)=-s;
end
pv=pv/sum(pv);
    
