function [V,l,lambdas,e,tail]=steroid(A,varargin)
% [V,l,lambdas,e,tail]=steroid(A,method) or [V,tail]=steroid(A,method)
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
% lambdas   =   vector, contains the weights in the STEROID,
%
% e         =   scalar, residual that is not described by the span of V,
%
% tail      =   tensor, symmetric tensor built up from the cross-product
%               contributions in the STEROID,
%
% A         =   tensor, symmetric d-way tensor,
%
% method    =   string, optional: determines how the least-squares problem 
%               W*d=vec(A) will be solved. Possible choices are:
%
%               'bigW': constructs the W matrix with n^d rows, only
%               feasible for small n and d. This is the default option,
%
%               'WtW': solves the much smaller W^T*W*d=W^T*vec(A) problem
%               at the cost of a worse condition number,
%
%               'Wsym': solves the smaller S*W*d=S*vec(A), where S is a row
%               selection matrix that selects only the nchoosek(d+n-1,n-1)
%               distinct entries in A.
%
% Reference
% ---------
%
% Symmetric Tensor Decomposition by an Iterative Eigendecomposition
% Algorithm
%
% 2014, 2015, Kim Batselier & Ngai Wong

n=size(A,1);
d=length(size(A));
doriginal=d;
a=A(:);

if sum(n(1)==size(A)) ~= d
	error('A needs to be a cubical tensor.');
end

numberoflevels=ceil(log2(d));
r=zeros(1,numberoflevels);
dtemp=d;
for i=1:length(r)
    if mod(dtemp,2)==0
        dtemp=dtemp/2;        
    else
        dtemp=(dtemp+1)/2;
    end
    r(i)=n^(dtemp); % number of branches per cluster at level i
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
embedtime=[];

if mod(d,2)==1
    % odd-order, need to embed
    tic
    A=embed(A);
    test=toc;
    embedtime=[embedtime test];
    d=(d+1)/2;
else
    d=d/2;
end

tic
% first eig
[V1,D1]=eig(reshape(A,[n^d n^d]));
test=toc;
eigtime=[eigtime test];
Dt{1}=diag(D1);
Vt{1}=V1;
% [Dt{1} I]=sort(abs(diag(D1)),'descend');
% Vt{1}=V1(:,I);
L{1}=diag(D1);
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are computing eigdecomp

if length(r)==1
    V=V1;   
else
    V = spalloc(n,prod(r),n*prod(r)); 
%     V=zeros(n,prod([sum(abs(Dt{1})>length(Dt{1})*eps(max(Dt{1}))) r(2:end)]));
    vcolcounter=1;
end
tol=n^d*eps(max(abs(Dt{whichvcounter})));
lambdas=kron(L{1}, ones(nleaf/length(L{1}),1));

for i=1:length(r)-1           % outer loop over the levels
    Llevel=[];
% 	tol=n^(d/(2^i))*eps(max(abs(Dt{whichvcounter})));    
    for j=1:prod(r(1:i))      % inner loop over the number of eigs for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if ~isempty(Dt{whichvcounter}) && abs(Dt{whichvcounter}(col)) > tol
            if mod(d,2)==1
                % odd-order, need to embed
                tic
                tempV=embed(reshape(Vt{whichvcounter}(:,col),n*ones(1,d)));
                test=toc;
                embedtime=[embedtime test];                dtemp=(d+1)/2;
            else
                tempV=reshape(Vt{whichvcounter}(:,col),n*ones(1,d));
                dtemp=d/2;
            end
            tic
            [V1,D1]=eig(reshape(tempV,[n^dtemp n^dtemp]));
            test=toc;
            eigtime=[eigtime test];

            Vt{counter}=V1;
            Dt{counter}=diag(D1);
            L{counter}=diag(D1).^(2^(i));
            Llevel=[Llevel;L{counter}];
            if i==length(r)-1                
                V(:,vcolcounter:vcolcounter+size(V1,2)-1)=V1;
                vcolcounter=vcolcounter+size(V1,2);
            end
        else
            L{counter}=zeros(n^dtemp,1);
            Llevel=[Llevel;L{counter}];
            if i==length(r)-1                
                vcolcounter=vcolcounter+size(V1,2);
            end            
        end
        counter=counter+1;
        if rem(j,length(Dt{whichvcounter}))==0
%             V{whichvcounter}=[];
            whichvcounter =  whichvcounter+1;
        end
    end
    d=dtemp;
    Llevel=kron(Llevel, ones(nleaf/length(Llevel),1));
    lambdas=lambdas.*Llevel;
end

clear A D1 V1 Dt Vt L Llevel col colcounter whichvcounter vcolcounter

% remove zero lambdas and corresponding vectors
V=full(V(:,abs(lambdas)>tol));
lambdas=lambdas(abs(lambdas)>tol);

if isempty(varargin)
    method='bigW';
else
    method=varargin{1};
end

if nargout==2
    % only V vectors and tail are required
    % need to compute the tail
    switch lower(method)
        case {'bigw'}
             W=zeros(n^doriginal,size(V,2));
             for i=1:size(W,2)
                 W(:,i)=mkron(V(:,i),doriginal);
             end
             l=reshape(a-W*lambdas,n*ones(1,doriginal));  % compute symmetric tail
        case {'wtw','wsym'}
            if length(lambdas) ~=size(V,2)
                disp('Warning: length lambdas not equal to number of V vectors');
                l=lambdas;
            else
                % compute symmetric tail
                head=zeros(n^doriginal,1);
                for i=1:length(lambdas)
                    head=head+lambdas(i)*mkron(V(:,i),doriginal);
                end
                l=reshape(a-head,n*ones(1,doriginal));
            end
    end
    return
end
  
            
%% solve the linear system W*d=vec(A)
switch lower(method)
    case {'bigw'}
        % original LS problem, no symmetry exploited
        W=zeros(n^doriginal,size(V,2));
        for i=1:size(W,2)
            W(:,i)=mkron(V(:,i),doriginal);
        end
        tic
        l=W\a;
        lstime=toc;
%         time(2)=toc;
        I=find(l);
        e=norm(a-W(:,I)*l(I));                          % compute residual
        tail=reshape(a-W*lambdas,n*ones(1,doriginal));  % compute symmetric tail
    case {'wtw',}
        % W^T*W
        WtW=(V'*V).^doriginal;
        % update righ-hand-side of LS problem, X^T*vec(A)
        b=zeros(size(V,2),1);
        for i=1:size(V,2)
            b(i,1)=mkron(V(:,i),doriginal)'*a;
        end
        tic
        l=WtW\b;
        lstime=toc;
%         time(2)=toc;
%         % solve linear system with SVD
%         [Uw Sw Vw]=svd(WtW);
%         s=diag(Sw);
%         rankW=sum(s>size(WtW,2)*eps(s(1)));        
%         Uw=Uw(:,1:rankW);
%         Sw=Sw(1:rankW,1:rankW);
%         Vw=Vw(:,1:rankW);
%         d=Vw*diag(1./diag(Sw))*Uw'*b;
        
        % compute residual
        ahat=zeros(n^doriginal,1);
        I=find(l);
        for i=1:length(I)
            ahat=ahat+l(I(i))*mkron(V(:,I(i)),doriginal);
        end        
        e=norm(a-ahat);
        
        % compute symmetric tail
        head=zeros(n^doriginal,1);
        for i=1:length(lambdas)
            head=head+lambdas(i)*mkron(V(:,i),doriginal);
        end   
        tail=reshape(a-head,n*ones(1,doriginal));
    case {'wsym'}
        % original Ls problem, symmetry exploited
        mons=getMonBase(doriginal,n);
        lindex=exp2ind(mons);
        b=a(lindex);
        Wsym=zeros(size(mons,1),size(V,2));
        for i=1:size(Wsym,1)
           Wsym(i,:)=prod(V.^(mons(i,:)'*ones(1,size(V,2))),1);
        end
        tic
        l=Wsym\b;
        lstime=toc;
%         time(2)=toc;
        I=find(l);
        e=norm(b-Wsym(:,I)*l(I));                          % compute residual
        
        % compute symmetric tail
        head=zeros(n^doriginal,1);
        for i=1:length(lambdas)
            head=head+lambdas(i)*mkron(V(:,i),doriginal);
        end   
        tail=reshape(a-head,n*ones(1,doriginal));
end
% disp(['error: ' num2str(e) ', number of eigs: ' num2str(length(eigtime)) ', number of embeddings: ' num2str(length(embedtime)) ', total eigtime: ' num2str(sum(eigtime)) ', total embedtime: ' num2str(sum(embedtime)) ',V time: ' num2str(sum(eigtime)+sum(embedtime)) ', LS time: ' num2str(lstime) ]);
% disp(['error: ' num2str(e) ', LS time: ' num2str(lstime) ]);
% disp(['error: ' num2str(e) ', V vectors: ' num2str(time(1)) ', LS: ' num2str(time(2)) ', total time: ' num2str(sum(time))])

end