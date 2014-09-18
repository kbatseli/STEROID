function e=symcheck(A)
% e=symcheck(A)
% -------------
% Checks whether a given tensor is symmetric.
%
% e         =   vector, e(i) contains ||A-A_sigma(i)||_2 with sigma(i)
% 				the i-th permutation of indices,
%
% A			=   tensor, arbitrary d-way tensor.

% Reference
% ---------
%
% 2014, Kim Batselier
d=length(size(A));
indices=perms([1:d]);

for i=1:size(indices,1)
   e(i)=norm(reshape(A-permute(A,indices(i,:)),[1 numel(A)])); 
end

end
