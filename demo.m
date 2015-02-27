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

