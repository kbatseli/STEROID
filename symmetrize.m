function As=symmetrize(A)
% As=symmetrize(A)
% ----------------
% Turns the tensor A into a symmetric tensor As.
%
% As        =   tensor, symmetric tensor obtained from A,
%
% A			=   tensor, arbitrary d-way tensor.

% Reference
% ---------
%
% 2014, Kim Batselier

d=length(size(A));
indices=perms([1:d]);
As=zeros(size(A));
for i=1:size(indices,1)
	As=As+permute(A,indices(i,:));
end
As=As/size(indices,1);

end
