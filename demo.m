% create a random 8-way symmetric tensor
A=randsymten(8,3);

% check whether it is symmetric
norm(symcheck(A))

% compute its 5-way approximation
[V5,D5,lambdas5,e5,X5,tail5]=steroid(A,5);

% compute its 8-way approximation
[V8,D8,lambdas8,e8,X8,tail8]=steroid(A);
