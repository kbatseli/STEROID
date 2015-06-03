% create a random 7-way symmetric tensor
clear all
A=randsymten(7,3);

% check whether it is symmetric
norm(symcheck(A))

% compute its STEROID with the original method
tic,
[V1,d1,e1]=steroid(A);
toc

% compute its STEROID with symmetry exploitation
tic,
[V2,d2,e2]=steroid(A,'wsym');
toc

% compute its STEROID with X^T*X
tic,
[V3,d3,e3]=steroid(A,'wtw');
toc

% compute possible Z-eigenpairs of A using PQRST
[lambda,V,err,itr]=pqrst(A);

% Suppose we have the following homogeneous polynomial
polyA{1,1}=[1 -2];
polyA{1,2}=[3 0 0;0 1 2];
% which corresponds to the following 3rd-order tensor of dimension 3
% A(1,1,1)=1;
% A(2,3,3)=A(3,2,3)=A(3,3,2)=-2;

% compute its STEROID for increasing R
[Vs,ds,es]=steroids(polyA,1);
% check residual 
es
[Vs,ds,es]=steroids(polyA,3);
es
[Vs,ds,es]=steroids(polyA,5);
es

% compute possible Z-eigenpairs of polyA using PQRST
[lambda,V,err,itr]=pqrst(polyA,1e-16,100,1,5);
