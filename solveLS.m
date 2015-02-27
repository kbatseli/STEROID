function [d e]=solveLS(V,A,varargin)
% [d e]=solveLS(V,A,method)
% -------------------------
% Tries to solve the least-squares problem W*d=vec(A), where each column of
% W is mkron(V(:,i),d).
%
% d         =   vector, contains the weights of each of the terms defined
%               by the columns of V in the decomposition,
%
% e         =   scalar, residual that is not described by the span of V,
%
% V         =   matrix, each column corresponds with a vector that
%               determines 1 rank-1 symmetric tensor,
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
doriginal=length(size(A));
a=A(:);

if isempty(varargin)
    method='bigW';
else
    method=varargin{1};
end

%% solve the linear system W*d=vec(A)
switch lower(method)
    case {'bigw'}
        % original LS problem, no symmetry exploited
        W=zeros(n^doriginal,size(V,2));
        for i=1:size(W,2)
            W(:,i)=mkron(V(:,i),doriginal);
        end
        d=W\a;
        I=find(d);
        e=norm(a-W(:,I)*d(I));                          % compute residual
    case {'wtw',}
        % W^T*W
        WtW=(V'*V).^doriginal;
        % update righ-hand-side of LS problem, X^T*vec(A)
        b=zeros(size(V,2),1);
        for i=1:size(V,2)
            b(i,1)=mkron(V(:,i),doriginal)'*a;
        end
        d=WtW\b;
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
        I=find(d);
        for i=1:length(I)
            ahat=ahat+d(I(i))*mkron(V(:,I(i)),doriginal);
        end        
        e=norm(a-ahat);              
    case {'wsym'}
        % original Ls problem, symmetry exploited
        mons=getMonBase(doriginal,n);
        lindex=exp2ind(mons);
        b=a(lindex);
        Wsym=zeros(size(mons,1),size(V,2));
        for i=1:size(Wsym,1)
           Wsym(i,:)=prod(V.^(mons(i,:)'*ones(1,size(V,2))),1);
        end
        d=Wsym\b;
        I=find(d);
        e=norm(b-Wsym(:,I)*d(I));                          % compute residual        
end


end