%%%
%Takes a matrix containing a variaty of values in each column vector
%and returns a matrix containing every combination between the elements
%in those vectors.
%%%

function grid=createGrid(gridValues,s)

k=size(gridValues,2);
temp=vectorCombinations(s);

grid=zeros(prod(s),k);
for i=1:k
	grid(:,i)=gridValues([temp(:,i)],i);
end

