function [U,V,D,e,X]=svdsteroid(A,ddesired)
% [U,V,D,e,X]=svdsteroid(A,desiredd)
% ---------------------------------------------
% Combination of an SVD and STEROID for a tensor A that is symmetric in all
% but the first mode. In addition, since STEROID is used, we further need that the
% order of the tensor minus one is a power of 2.
%
% U         =   matrix, each column corresponds with a vector that
%               determines the first mode of A in the reconstruction,
%
% V         =   matrix, each column corresponds with a vector that
%               determines 1 rank-1 symmetric tensor for all modes except the first,
%
% D         =   vector, contains the weights of each of the terms defined
%               by the columns of U,V in the decomposition,
%
% e         =   scalar, residual that is not described by the span of U,V,
%
% X         =   matrix, each column contains the Kronecker product of the
%               corresponding vectors in U,V,
%
% A         =   tensor,  symmetric in all but the first mode,
%
% desiredd  =   scalar, optional: in case A is the embedding of a symmetric cubical
%               tensor of lower order, then providing the original order
%               will return the decomposition of the original tensor.
%               Default value = the order of A.               
%
% Reference
% ---------
%
% 2015, Kim Batselier

% check whether d-1 is a power of 2
n=size(A);
d=length(n);
if ceil(log2(d-1))-log2(d-1) ~=0
	error('The order of A needs to be a power of 2.');
end

if nargin==2
	ddesired=varargin{1};
else
	ddesired=d;
end

% first step is svd
[U S1 V1]=svd(reshape(A,[n(1) prod(n(2:end))]),'econ');
X=[];
for i=1:size(V1,2)
    [V{i},Ds{i},lambdas{i},es{i},Xs{i},tails{i}]=steroid(reshape(V1(:,i),n(2:end)),ddesired-1);
    X=[X mkron(Xs{i},U(:,i))];
end

% retrieve the "eigenvalues"
if ddesired < d
    cmpstring='A=A(:';
    for i=2:ddesired
        cmpstring=[cmpstring ',:'];
    end
    for i=ddesired+1:d
        cmpstring=[cmpstring ',1'];
    end
    cmpstring=[cmpstring ');'];
    eval(cmpstring);
end
D=X\reshape(A,[n(2)^ddesired 1]);
e=norm(reshape(A,[n(2)^ddesired 1])-X*D);


end