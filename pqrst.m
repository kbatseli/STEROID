function [lambda,V,err,itr]=pqrst(A,varargin)
%  [lambda,V,err,itr]=pqrst(A,tol,kmax,shifting,R)
% ------------------------------------------------
% Tries to compute Z-eigenpairs of the symmetric tensor A using the QR
% algorithm on matrix slices of permuted versions of A.
%
% lambda    =   vector, each entry  corresponds with a possible
%               Z-eigenvalue,
%
% V         =   matrix, each column corresponds with a possible Z-eigenvector,
%
% err       =   vector, contains ,
%
% e         =   scalar, residual that is not described by the span of V,
%
% tail      =   tensor, symmetric tensor built up from the cross-product
%               contributions in the STEROID,
%
% A         =   tensor, symmetric d-way tensor,
%
%
% Reference
% ---------
%
% A QR algorithm for symmetric tensors
%
% 2014, 2015, Kim Batselier & Ngai Wong

switch length(varargin)
    case 0
        % default tolerance
        tol=1e-10;
        % default maximal number of iterations
        kmax=100;
        % default shift
        shift =0;
    case 1
        tol=varargin{1};
        % default maximal number of iterations
        kmax=100;
        % default shift
        shift =0;
    case 2
        tol=varargin{1};
        % default maximal number of iterations
        kmax=varargin{2};
        % default shift
        shift =1;        
    case 3
        tol=varargin{1};
        % default maximal number of iterations
        kmax=varargin{2};
        % default shift
        shift =varargin{3};
    case 4
        tol=varargin{1};
        % default maximal number of iterations
        kmax=varargin{2};
        % default shift
        shift =varargin{3};
        % number of eigenvectors to compute for each eigs call in steroids
		R = varargin{4};
    otherwise
        error('No more than 4 additional input arguments are supported.');
end

if ~iscell(A)
	% A is a d-way symmetric tensor
	n=size(A); % gets the dimension n for each way
	d=length(n);
	n=n(1);
	[V1,l1,e1]=steroid(A);
	if e1>n^d*eps(norm(A(:)))
	    error('Steroid is not exact. Need an additional iteration on the tail tensor.');
	end
else
	% A is sparse symmetric tensor given as a homogeneous polynomial
	% check whether R argument was given
	if length(varargin)~=4
		error('Please provide the 4th optional R argument.')
	end
	n=size(A{1,2},2);
	d=sum(A{1,2}(1,:));
	[V1,l1,e1]=steroids(A,R);
	if e1>n^d*eps(norm(A{1,2}))
	    error('Steroid is not exact. Need an additional iteration on the tail tensor or larger R.');
	end
end

p=perms(n:-1:1);
E=eye(n);

lambda=zeros(1,n*length(p)); % holds eigenvalues
V=zeros(n,n*length(p)); % holds eigenvectors
err=zeros(1,n*length(p)); % holds error in QR algo
itr=zeros(1,n*length(p)); % holds iteration count for each slice

for i=1:size(p,1)
    U=E(:,p(i,:)); %similarity transform matrix holder
    for j=1:n
        ej=zeros(n,1);
        ej(j)=1;
        for k=1:kmax
            M=U'*V1*diag(l1'.*(U(:,j)'*V1).^(d-2))*V1'*U;
            epsilon=norm(M(:,j)-M(j,j)*ej);
            if epsilon<tol
                break;
            end
            s=0;
            if shift
                lambdamin=min(real(eig(M)));
                if lambdamin <0
                    s=-lambdamin+shift;
                end
            end
            [Q,R]=qr(M+s*eye(n,n));
            U=U*Q;
        end
        
        lambda((i-1)*n+j)=M(j,j);
        V(:,(i-1)*n+j)=U(:,j);
        err((i-1)*n+j)=norm( V1*((V1'*U(:,j)).^(d-1).*l1)-U(:,j)*M(j,j) );
        itr((i-1)*n+j)=k;
    end
end

end
