function [V,D,lambdas,e,X,tail]=steroid(A,varargin)
% [V,D,lambdas,e,X,tail]=steroid(A,desiredd)
% ------------------------------------------
% Symmetric Tensor Eigen Rank-One Iterative Decomposition. Decomposes a
% symmetric cubical tensor into a real linear combination of real rank-1 symmetric tensors.
%
% V         =   matrix, each column corresponds with a vector that
%               determines 1 rank-1 symmetric tensor,
%
% D         =   vector, contains the weights of each of the terms defined
%               by the columns of V in the decomposition,
%
% lambdas   =   vector, contains the weights in the STEROID,
%
% e         =   scalar, residual that is not described by the span of V,
%
% X         =   matrix, each column contains the Kronecker product of the
%               corresponding vector in V,
%
% A         =   tensor, cubical symmetric tensor,
%
% tail      =   tensor, symmetric tensor built up from the cross-product
%               contributions in the STEROID,
%
% desiredd  =   scalar, optional: in case A is the embedding of a symmetric cubical
%               tensor of lower order, then providing the original order
%               will return the decomposition of the original tensor.
%               Default value = the order of A.               
%
% Reference
% ---------
%
% 2014, Kim Batselier

n=size(A,1);
d=length(size(A));
if ceil(log2(d))-log2(d) ~=0
	error('The order of A needs to be a power of 2.');
else
	numberoflevels=log2(d);
end

if nargin==2
	ddesired=varargin{1};
else
	ddesired=d;
end

if sum(n(1)==size(A)) ~= d
	error('A needs to be a cubical tensor.');
end		

r=n.^(d./2.^[1:numberoflevels]);
eigsperlevel=ones(1,numberoflevels);
totaleigs=0;
for i=1:length(r)-1
    eigsperlevel(i+1)=prod(r(1:i));
    totaleigs=totaleigs+eigsperlevel(i+1); 
end
nleaf=prod(r);

Vt=cell(1,totaleigs);
Dt=cell(1,totaleigs);
L=cell(1,totaleigs);

% first eig
[V1,D1]=eig(reshape(A,[n^(d/2) n^(d/2)]));
Vt{1}=V1;
Dt{1}=diag(D1);
L{1}=diag(D1);
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are computing eigdecomp
V=[];
X=[];
for i=1:length(r)-1           % outer loop over the levels
	tol=n^d*eps(max(abs(Dt{whichvcounter})));
% 	tol=n^(d/(2^i))*eps(max(abs(Dt{whichvcounter})));    
    for j=1:prod(r(1:i))      % inner loop over the number of eigs for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if ~isempty(Dt{whichvcounter}) && abs(Dt{whichvcounter}(col)) > tol
                [V1,D1]=eig(reshape(Vt{whichvcounter}(:,col),[n^(d/2^(i+1)) n^(d/2^(i+1))]));
                Vt{counter}=V1;
                Dt{counter}=diag(D1);
                L{counter}=diag(D1).^(2^(i));
                if i==length(r)-1
                    V=[V V1];
                    for k=1:n
                        X=[X mkron(V1(:,k),ddesired)];
                    end
                end
        else
           L{counter}=zeros(n^(d/2^(i+1)),1); 
        end
        counter=counter+1;
        if rem(j,length(Dt{whichvcounter}))==0
%             V{whichvcounter}=[];
            whichvcounter =  whichvcounter+1;
        end
    end
%     whichvcounter = whichvcounter+1;
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
D=X\reshape(A,[n^ddesired 1]);
e=norm(reshape(A,[n^ddesired 1])-X*D);

% compute the lambdas
Llevel=cell(1,length(r));   % cat each level singular values into 1 vector
counter=1;
for i=1:length(r),
    for j=1:eigsperlevel(i),
        Llevel{i}=[Llevel{i}; L{counter}];
        counter=counter+1;
    end
end

for i=1:length(r),             % make all singular value vectors the same size (number of leaves)
    Llevel{i}=kron(Llevel{i}, ones(nleaf/length(Llevel{i}),1));
end

lambdas=ones(nleaf,1);         % output singular values at each leaf
for i=1:length(r),
    lambdas=lambdas.*Llevel{i};
end

% compute the symmetric tail part
head=reshape(X*lambdas(find(lambdas)),size(A));
tail=A-head;

end
