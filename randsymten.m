function A=randsymten(d,n)
% A=randsymten(d,n)
% -----------------
% Generates a random d-way cubical symmetric tensor of dimension n.
%
% A         =   tensor, random symmetric cubical tensor,
%
% d			=   scalar, number of modes (order),
%
% n 		=	scalar, dimension of each of the modes.
%
% Reference
% ---------
%
% 2014, Kim Batselier

A0=randn(n*ones(1,d));
A=zeros(size(A0));
indices=perms([1:d]);
for i=1:size(indices,1)
	A=A+permute(A0,indices(i,:));
end
A=A/size(indices,1);	

end
