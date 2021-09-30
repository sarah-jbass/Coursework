function [b,o,j,conv]=repeatSearch(func,b,options,data,N,convTOL);

oDiff=convTOL+20;
oldo=0;
i=1;
conv=0;

while i<=N & oDiff>convTOL
   [b,o,j]=fminsearch(func,b,options,data);
   oDiff=abs(o-oldo);
   oldo=o;
   i=i+1;
end

if oDiff<=convTOL
   conv=1;
end



