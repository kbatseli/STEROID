function A=poly2mat(polyA)
% A=poly2mat(polyA)
% -----------------
% Converts a homogeneous polynomial (in polysys cell format) of even degree
% d into a square n^(d/2) sparse matrix.
%
% A         =   sparse matrix, symmetric tensor reshaped into n^(d/2)
%               square matrix,
%
% polyA     =   cell, polysys cell for 1 homogeneous polynomial,
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
if sum(alld==max(alld)) == length(alld) && mod(alld(1),2)==0
    d=alld(1);
else
    error('polyA needs to be a homogeneous polynomial and of even order.')
end

rowI=[];
colI=[];
val=[];
nzmax=0;

for i=1:length(polyA{1,1})
    [lindex index]=exp2ind(polyA{1,2}(i,:));
    % check how many entries we need for this coefficient
    indices=intersect(perms(index),perms(index),'rows');
    numberofentries = size(indices,1);
    nzmax=nzmax+numberofentries;
       
    % collect all values
    val=[val;polyA{1,1}(i)*ones(numberofentries,1)];
    
    % compute corresponding row and column indices
    for j=1:numberofentries
        rowindex=indices(j,1:d/2);
        colindex=indices(j,d/2+1:end);
        rowtemp(j,1)=rowindex(1);
        coltemp(j,1)=colindex(1);
        for k=2:length(rowindex)
            rowtemp(j,1)=rowtemp(j,1)+(rowindex(k)-1)*n^(k-1);
            coltemp(j,1)=coltemp(j,1)+(colindex(k)-1)*n^(k-1);
        end        
    end
    rowI=[rowI;rowtemp];
    colI=[colI;coltemp];
    rowtemp=[];
    coltemp=[];
end
A=sparse(rowI,colI,val,n^(d/2),n^(d/2),nzmax);

end