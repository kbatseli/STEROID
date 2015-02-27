function B=embed(A)
%  B=embed(A)
% -----------
% Embeds a symmetric d-way tensor A where d is odd into a symmetric d+1-way
% tensor B such that B(:,:,...,:,1)=A.
%
% B         =   tensor, symmetric d+1-way tensor,
%
% A         =   tensor, symmetric d-way tensor with d odd.
%
%
% Reference
% ---------
%
% Symmetric Tensor Decomposition by an Iterative Eigendecomposition
% Algorithm
%
% 2015, Kim Batselier & Ngai Wong

n=size(A);
d=length(n);
n=n(1);

if mod(d,2) == 0
    disp('Given tensor A is already of even order.');
    B=[];
    return
end

[lindex indices]=exp2ind(getMonBase(d,n));
% extend indices to even order
indices=[indices ones(size(indices,1),1)];

B=zeros(prod(n*ones(1,d+1)),1);
   
for i=1:size(indices,1)
    temp=perms(indices(i,:));
    temp=intersect(temp,temp,'rows');
    for j=1:size(temp,1)
        string=['index=sub2ind(' num2str(n) '*ones(1,' num2str(d+1) ')'];
        for k=1:d+1
            string=[string ',temp(:,' num2str(k) ')'];
        end
        string=[string ');'];
        eval(string);
        B(index)=A(lindex(i));
    end
end

B=reshape(B,[n*ones(1,d+1)]);

end
