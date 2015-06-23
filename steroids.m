function [V,l,e,varargout]=steroids(polyA,R)
% [V,l,e,tail,lambdas]=steroids(polyA,R) or [V,tail]=steroids(polyA,R)
% --------------------------------------------------------------------
% Symmetric Tensor Eigen Rank-One Iterative Decomposition. Decomposes a
% symmetric tensor into a real linear combination of real rank-1 symmetric 
% tensors. If STEROID is called with only two output arguments, then only
% the V vectors and the symmetric tail will be computed and returned.
%
% V         =   matrix, each column corresponds with a vector that
%               determines 1 rank-1 symmetric tensor,
%
% l         =   vector, contains the weights of each of the terms defined
%               by the columns of V in the decomposition,
%
% e         =   scalar, residual that is not described by the span of V,
%
% tail      =   tensor, symmetric tensor built up from the cross-product
%               contributions in the STEROID,
%
% lambdas   =   vector, contains the weights in the STEROID,
%
% polyA     =   cell, polysys cell for 1 homogeneous polynomial,
%
% R         =   scalar, number of eigenvectors that need to be computed each eigs.
%
% Reference
% ---------
%
% Symmetric Tensor Decomposition by an Iterative Eigendecomposition
% Algorithm
%
% 2014, 2015, Kim Batselier & Ngai Wong

n=size(polyA{1,2},2);
alld=sum(polyA{1,2},2);
if sum(alld==max(alld)) == length(alld)
    d=alld(1);
else
    error('polyA needs to be a homogeneous polynomial.')
end
doriginal=d;
a=hpoly2vec(polyA);

numberoflevels=ceil(log2(d));
r=zeros(1,numberoflevels);
dtemp=d;
for i=1:length(r)
    if mod(dtemp,2)==0
        dtemp=dtemp/2;        
    else
        dtemp=(dtemp+1)/2;
    end
    if n^dtemp > R
        r(i)=R; % number of branches per cluster at level i
    else
        r(i)=n^dtemp;
    end
end

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

eigtime=[];

if mod(d,2)==1
    % odd-order, need to embed
    polyA{1,2}=polyA{1,2}+[ones(length(polyA{1,1}),1) zeros(length(polyA{1,1}),n-1)];
    d=(d+1)/2;
else
    d=d/2;
end

tic
% first eig
[V1,D1,flag]=eigs(poly2mat(polyA),R);
test=toc;
eigtime=[eigtime test];
if flag
    warning('Some eigenvalues did not converge.');
end
Dt{1}=diag(D1);
% set numerically zero values to zero and store the eigenvectors as sparse
% vectors again, we assume here that the eigenvectors are also sparse
tol=n^d*eps(max(abs(Dt{1})));
V1(abs(V1)<tol)=0;
Vt{1}=sparse(V1);
% [Dt{1} I]=sort(abs(diag(D1)),'descend');
% Vt{1}=V1(:,I);
L{1}=diag(D1);
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are computing eigdecomp

if length(r)==1
    V=V1;   
else
    V = spalloc(n,R^numberoflevels,n*R^numberoflevels); 
%     V=zeros(n,prod([sum(abs(Dt{1})>length(Dt{1})*eps(max(Dt{1}))) r(2:end)]));
    vcolcounter=1;
end

lambdas=kron(L{1}, ones(nleaf/length(L{1}),1));

for i=1:length(r)-1           % outer loop over the levels
    Llevel=[];
    for j=1:prod(r(1:i))      % inner loop over the number of eigs for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if ~isempty(Dt{whichvcounter}) && abs(Dt{whichvcounter}(col)) > tol
            % convert V vector into homogeneous polynomial
            polyV=vec2hpoly(Vt{whichvcounter}(:,col),d,n);
            if mod(d,2)==1
                % odd-order, need to embed
                 polyV{1,2}=polyV{1,2}+[ones(length(polyV{1,1}),1) zeros(length(polyV{1,1}),n-1)];
                 dtemp=(d+1)/2;
            else
                dtemp=d/2;
            end
            tempV=poly2mat(polyV);
            tic
            [V1,D1]=eigs(tempV,min(size(tempV,1),R));
            test=toc;
            eigtime=[eigtime test];
            if flag
                warning('Some eigenvalues did not converge.');
            end
            V1(abs(V1)<tol)=0;
            Vt{counter}=sparse(V1);;
            Dt{counter}=diag(D1);
            L{counter}=diag(D1).^(2^(i));
            Llevel=[Llevel;L{counter}];
            if i==length(r)-1                
                V(:,vcolcounter:vcolcounter+size(V1,2)-1)=V1;
                vcolcounter=vcolcounter+size(V1,2);
            end
        else
            L{counter}=zeros(min(size(tempV,1),R),1);
            Llevel=[Llevel;L{counter}];
            if i==length(r)-1                
                vcolcounter=vcolcounter+size(V1,2);
            end            
        end
        counter=counter+1;
        if rem(j,n^d)==0
            whichvcounter =  whichvcounter+1;
        end
    end
    d=dtemp;
    Llevel=kron(Llevel, ones(nleaf/length(Llevel),1));
    lambdas=lambdas.*Llevel;
end

clear D1 V1 Dt Vt L Llevel col colcounter whichvcounter vcolcounter

% remove zero lambdas and corresponding vectors
V=full(V(:,abs(lambdas)>tol));
lambdas=lambdas(abs(lambdas)>tol);

if nargout==2
    % only V vectors and tail are required
    % need to compute the tail
    head=zeros(n^doriginal,1);
    for i=1:length(lambdas)
        head=head+lambdas(i)*mkron(V(:,i),doriginal);
    end
    l=reshape(a-head,n*ones(1,doriginal));
    return
end
  
            
%% solve the linear system W*d=vec(A)
% Always use W^T*W
WtW=(V'*V).^doriginal;
% update righ-hand-side of LS problem, X^T*vec(A)
b=zeros(size(V,2),1);
for i=1:size(V,2)
    b(i,1)=mkron(V(:,i),doriginal)'*a;
end
tic
l=pinv(WtW)*b;
lstime=toc;

if nargout > 3
    % compute symmetric tail
    head=zeros(n^doriginal,1);
    for i=1:length(lambdas)
        head=head+lambdas(i)*mkron(V(:,i),doriginal);
    end
    varargout{1}=reshape(a-head,n*ones(1,doriginal));
    varargout{2}=lambdas;
else
    varargout=cell(1,1);
end
% compute residual
ahat=zeros(n^doriginal,1);
% remove zero lambdas and corresponding vectors
V=full(V(:,abs(l)>tol));
l=l(abs(l)>tol);
for i=1:length(l)
    ahat=ahat+l(i)*mkron(V(:,i),doriginal);
end        
e=norm(a-ahat);

% disp(['error: ' num2str(e) ', number of eigs: ' num2str(length(eigtime)) ', number of embeddings: ' num2str(length(embedtime)) ', total eigtime: ' num2str(sum(eigtime)) ', total embedtime: ' num2str(sum(embedtime)) ',V time: ' num2str(sum(eigtime)+sum(embedtime)) ', LS time: ' num2str(lstime) ]);
% disp(['error: ' num2str(e) ', LS time: ' num2str(lstime) ]);
% disp(['error: ' num2str(e) ', V vectors: ' num2str(time(1)) ', LS: ' num2str(time(2)) ', total time: ' num2str(sum(time))])

end
