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


if mod(d,2) == 0
    disp('Given tensor A is already of even order.');
    B=[];
    return
end

% first determine nonzero entries
I=find(A(:));
% convert linear indices to d-way indices
tempindices=zeros(size(I,1),d);
for i=1:size(I,1)
    string='[';
    for j=1:d-1
        string=[string 'tempindices(' num2str(i) ',' num2str(j) '),'];        
    end
    string=[string 'tempindices(' num2str(i) ',' num2str(d) ')]=ind2sub(n,I(' num2str(i) '));'];
    eval(string);
end
% now remove duplicate indices
indices=[];
lindex=[];
while ~isempty(tempindices)
    % store the index of a distinct entry
    indices=[indices;tempindices(1,:)];
    % keep the distinct entry value
    string=['A(' num2str(tempindices(1,1))];
    for i=2:d
        string=[string ',tempindices(1,' num2str(i) ')'];
    end
    string=[string ');'];   
    lindex=[lindex;eval(string)];
        
    temp=perms(tempindices(1,:));
    temp=intersect(temp,temp,'rows');
    tempindices=setdiff(tempindices,temp,'rows');
end
% extend indices to even order
indices=[indices ones(size(indices,1),1)];
n=n(1);
B=zeros(n^(d+1),1);    

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
        B(index)=lindex(i);
    end
end

B=reshape(B,[n*ones(1,d+1)]);

end
