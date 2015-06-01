% create a random 7-way symmetric tensor
A=randsymten(7,3);

% check whether it is symmetric
norm(symcheck(A))

% compute its STEROID with the original method
tic,
[V1,d1,lambdas1,e1,tail1]=steroid(A);
toc

% compute its STEROID with symmetry exploitation
tic,
[V2,d2,lambdas2,e2,tail2]=steroid(A,'wsym');
toc

% compute its STEROID with X^T*X
tic,
[V3,d3,lambdas3,e3,tail3]=steroid(A,'wtw');
toc

% Suppose we have the following homogeneous polynomial
polyA{1,1}=[1 -2];
polyA{1,2}=[3 0 0;0 1 2];
% which corresponds to the following 3rd-order tensor of dimension 3
% A(1,1,1)=1;
% A(2,3,3)=A(3,2,3)=A(3,3,2)=-2;

% compute its STEROID for increasing R
[Vs,ds,lambdass,es,tails]=steroid(polyA,1);
% check residual 
es
[Vs,ds,lambdass,es,tails]=steroid(polyA,3);
es
[Vs,ds,lambdass,es,tails]=steroid(polyA,5);
es
