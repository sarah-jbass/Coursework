function CovFn=R1_CovFn(GMat,nei_id,KMN)
% Calculate the variance-covariance matrix

N=size(GMat,1);
NG=size(GMat,2);
S=GMat; S=[S; zeros(1,NG)];
CovFn=zeros(NG,NG); 
for i=1:N
    t1=GMat(i,:);
    t2=nei_id(i,:);
    t3= sum(S(t2,:))*KMN + t1;
    CovFn = CovFn + t1'*t3;
end;
CovFn=CovFn/N;
