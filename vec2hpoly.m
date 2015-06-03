function polyV=vec2hpoly(V,d,n)
% polyV=vec2hpoly(V,d,n)
% ----------------------
% Converts a sparse vector into a homogeneous polynomial (in polysys cell format).
%
% polyV     =   cell, polysys cell for 1 homogeneous polynomial,
%
% V         =   sparse vector, coefficient vector of a homogeneous
%               polynomial,
%
% d         =   scalar, degree of the homogeneous polynomial,
%
% n         =   scalar, number of variables in the homogeneous polynomial.
%
% Reference
% ---------
%
% Symmetric Tensor Decomposition by an Iterative Eigendecomposition
% Algorithm
%
% 2014, 2015, Kim Batselier & Ngai Wong
polyV=cell(1,2);
I=find(V);
tempindices=zeros(size(I,1),d);
for i=1:length(I)
    string='[';
    for j=1:d-1
        string=[string 'tempindices(' num2str(i) ',' num2str(j) '),'];        
    end
    string=[string 'tempindices(' num2str(i) ',' num2str(d) ')]=ind2sub(n*ones(1,d),I(' num2str(i) '));'];
    eval(string);
end
% now remove duplicate indices
indices=[];
lindex=[];
while ~isempty(tempindices)
    % store the index of a distinct entry    
    polyV{1,1}=[polyV{1,1},V(I(1))];
    % store the exponent of that entry
    exponent=zeros(1,n);
    for i=1:n
        exponent(i)=sum(tempindices(1,:)==i);
    end
    polyV{1,2}=[polyV{1,2};exponent];
    % keep the distinct entry value
           
    temp=perms(tempindices(1,:));
    temp=intersect(temp,temp,'rows');
    [tempindices,IA]=setdiff(tempindices,temp,'rows');
    allI=1:length(I);
    allI(IA)=[];
    V(I(allI))=0;
    I(allI)=[];
end
end