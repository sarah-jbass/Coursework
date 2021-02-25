function [Ktraj, Ctraj, diff] = trajectory(kss,css,D,palpha,psigma,pbeta,pdelta)
C0 = 3.82:0.0000001:css;
C = C0;
K = kss+0*C;
nd = length(D);
C = repmat(C,nd,1); %takes starting point and turns them into a matrix with 700 rows (months)
K = repmat(K,nd,1);
for i=1:nd-1
    [K(i+1,:),C(i+1,:)] = iterate(K(i,:),C(i,:),D(i+1),psigma,palpha,pbeta,pdelta);
end
[diff,ind] = min(abs(C(end,:)-css)+abs(K(end,:)-kss));
Ctraj = C(:,ind)';
Ktraj = K(:,ind)';
end
