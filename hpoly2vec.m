function a=hpoly2vec(polyA)
% a=hpoly2vec(polyA)
% ------------------
% Converts a homogeneous polynomial (in polysys cell format) into the
% vectorization of the corresponding symmetric tensor.
%
% a         =   sparse vector, vectorization of symmetric tensor
%               corresponding with the homogeneous polynomial,
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
if sum(alld==max(alld)) == length(alld)
    d=alld(1);
else
    error('polyA needs to be a homogeneous polynomial.')
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
    val=[val;polyA{1,1}(i)*ones(numberofentries,1)/numberofentries];
    colI=[colI;ones(numberofentries,1)];
    
    % compute corresponding row and column indices
    for j=1:numberofentries
       string=['rowI=[rowI;sub2ind(n*ones(1,d),indices(j,1)'];
       for k=2:d
           string=[string ',indices(j,' num2str(k) ')' ];
       end
       eval([string ')];' ])
    end
end
a=sparse(rowI,colI,val,n^d,1,nzmax);

end
