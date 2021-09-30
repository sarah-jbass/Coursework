%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes n the rows of a matrix, and k the columns of a matrix
%and creates a matrix with all permutations of combinations
%between the vectors in each column.  Hence, it creates a matrix
%of size n^k X k.

function f=vectorCombinations(s)

gridSize=prod(s);
rowsS=size(s,1);
f=zeros(gridSize,rowsS);

for i=1:rowsS
   f(:,i)=kron(ones(prod(s(i:rowsS))/s(i),1),kron([1:s(i)]',ones(prod(s(1:i))/s(i),1)));   
end