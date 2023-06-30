function matrixset = matrixgenerator(n,p,seed,classes,numbers)
%%MATRIXGENERATOR produces sets of stochastic matrices of which we want to
% compute a pth root.
% INPUT: n size of the sample matrices
%        seed some of the matrices are generated by using a random number
%        generator, you can select here the seed to the generator
%        classes: array of integers specifying what class of matrix has to
%        be generated; see the following for a description of the classes.
%        numbers: how many arrays per class we need to generate, the vector
%           must have the same length as the vector of the classes.
% OUTPUT: cell-array containing the generated matrices

% SANITY CHECKS
if n < 0
    error("MATRIXGENERATOR:negative size n")
end
if ~exist('seed', 'var') || isempty(seed)
    seed = 42;
end
if length(classes) ~= length(numbers)
    error("MATRIXGENERATOR:classes and number of different length")
end

matrixset = cell(sum(numbers),1);
% Fix the seed
rng(seed);
% Class names
names = {'Uniform random',sprintf('%dth power of Uniform random',p),...
    'Exponentials of intensity matrix','K80 (embeddable)',...
    'K80 (not embeddable)','Pei Matrix (embeddable)'};
classessizes = [n,n,n,4,4,n];

matrixcounter = 1;
for i=1:length(classes)
    for j=1:numbers(i)
        switch classes(i)
            case 1
                % Random n x n matrices with elements from the uniform
                % distribution on [0,1] which are then scaled to a
                % stochastic matrix by dividing each element by its
                % row sum.
                A = rand(n,n);
                D = diag(sum(A,2));
                matrixset{matrixcounter} = D\A;
            case 2
                % X^p matrices where X is an n x n matrices with elements
                % from the uniform distribution on [0,1] which are then
                % scaled to a stochastic matrix by dividing each element
                % by its row sum and p given in input
                A = rand(n,n);
                D = diag(sum(A,2));
                matrixset{matrixcounter} = mpower(D\A,p);
            case 3
                % A = expm(Q) where Q is an intensity matrix obtained by 
                % generating a random n x n matrix with elements from the 
                % uniform distribution on [0, 1] and then adjusting the 
                % diagonal elements such that each row sum is zero.
                A = rand(n,n);
                A(1:1+n:end) = 0;
                D = diag(sum(A,2));
                matrixset{matrixcounter} = expm(A - D);
            case 4
                % K80 Model (embeddable). With this choice of parameters
                % the resulting matrix is always embeddable
                b = rand(1);
                c = sqrt(b) - b;
                a = 1 -b-2*c;
                A0 = [a,b;a,b];
                E = ones(2,2);
                matrixset{matrixcounter} = [A0,c*E;c*E,A0];
            case 5
                % K80 Model (not embeddable). With this choice of
                % parameters the matrix is always not embeddable
                b = rand(1);
                c = (1-2*b)/2;
                a = 1 -b-2*c;
                A0 = [a,b;a,b];
                E = ones(2,2);
                matrixset{matrixcounter} = [A0,c*E;c*E,A0];
            case 6
                % Pei Matrix
                % Section 2.6.4 of Lin thesis for this choice of parameters
                % they always have a primary stochastic pth root.
                I = eye(n,n);
                J = ones(n,n);
                alpha = rand(1)-(1/(n-1))^p;
                beta = (1 - alpha)/n;
                matrixset{matrixcounter} = alpha*I+beta*J;
            otherwise
                error("MATRIXGENERATOR:Uknown class");
        end
        matrixcounter = matrixcounter + 1;
    end
end

fprintf("Completed generation of %d matrices (seed = %d):\n", ...
    matrixcounter-1,seed);
for i=1:length(classes)
    fprintf("\t- %d matrices of class %s (size %d x %d)\n", ...
        numbers(i),names{classes(i)},classessizes(classes(i)),...
        classessizes(classes(i)));
end



end