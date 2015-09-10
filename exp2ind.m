function [lindex index]=exp2ind(exponent)
% [lindex index]=exp2ind(exponent)
% --------------------------------
% Converts homogeneous monomial exponents into a linear and tensor index of a d-way
% symmetric tensor of dimension n.
%
% lindex    =   vector, each entry is linear index corresponding with
%               monomial obtained from exponent(k,:),
%
% index 	=	matrix, each row index(k,:) is a d-way index of
%               corresponding linear index lindex(k),
%
% exponent  =   matrix, each row exponent(k,:) corresponds with an exponent
%               of a homogeneous monomial in n variables of degree d.
%
% Reference
% ---------
%
% 2015, Kim Batselier & Ngai Wong

[p n]=size(exponent);
d=sum(exponent(1,:));

index=zeros(p,d);

for i=1:p
	counter=1;
	for j=1:n
		if exponent(i,j)~=0
			for k=1:exponent(i,j)
				index(i,counter)=j;
				counter=counter+1;
			end
		end
	end
end
if d==1
	string=['lindex=sub2ind([' num2str(n) '*ones(1,' num2str(d) '),1]'];
else
	string=['lindex=sub2ind(' num2str(n) '*ones(1,' num2str(d) ')'];
end
for i=1:d
	string=[string ',index(:,' num2str(i) ')'];
end
string=[string ');'];
eval(string);
end
