function [U,V,d,e]=svdsteroid(A)
% [U,V,d,e]=svdsteroid(A)
% -----------------------
% Combination of an SVD and STEROID for a tensor A that is symmetric in all
% but the first mode.
%
% U         =   matrix, each column corresponds with a vector that
%               determines the first mode of A in the reconstruction,
%
% V         =   matrix, each column corresponds with a vector that
%               determines 1 rank-1 symmetric tensor for all modes except the first,
%
% d         =   vector, contains the weights of each of the terms defined
%               by the columns of U,V in the decomposition,
%
% e         =   scalar, residual that is not described by the span of U,V,
%
% A         =   tensor,  symmetric in all but the first mode,
%
% Reference
% ---------
%
% 2015, Kim Batselier

n=size(A);
d=length(n);

% first step is svd
[U S1 V1]=svd(reshape(A,[n(1) prod(n(2:end))]),'econ');
% now determine number of nonzero singular values
s1=diag(S1);
tol=prod(n(2:end))*eps(s1);
r=sum(s1>tol);
U=U(:,1:r);
X=[];
for i=1:r
    [V{i},tails{i}]=steroid(reshape(V1(:,i),n(2:end)));
    W=zeros(n(2)^(d-1),size(V{i},2));
    for j=1:size(V{i},2)
        W(:,j)=mkron(V{i}(:,j),d-1);
    end
    X=[X kron(W,U(:,i))];
end

% retrieve the "eigenvalues"
d=X\A(:);
e=norm(A(:)-X*d);


end
