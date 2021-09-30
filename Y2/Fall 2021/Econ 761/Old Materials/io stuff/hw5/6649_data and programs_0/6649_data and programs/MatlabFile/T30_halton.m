function Eps=T30_halton(prime,NR,NObs,N0)

k=prime;
s=[];
for j=1:k
    s=[s;(j-1)/k];
end;

for i=2:ceil(log(NR*NObs)/log(k))+1
    t=s;
    for j=1:k-1
        s=[s;t+j/k^i];
    end;
end;
T=zeros(NR,NObs);
T(:)=s(N0+1:NR*NObs+N0);
T=T'; Eps = norminv(T);

