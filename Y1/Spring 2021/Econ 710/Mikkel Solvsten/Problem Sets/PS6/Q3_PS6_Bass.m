clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\Y1\Spring 2021\Econ 710\Mikkel Solvsten\Problem Sets\PS6'

%Written by Sarah Bass, 3/3/21

%Params
beta0 = 1;
pdelta = [0 1 1 1]';
phis = [0 0.8];
ns = [40 70 100];
T = 4;
b_ols = zeros(length(ns),length(phis),5);
b_fe = zeros(length(ns),length(phis),4);
CIr = zeros(length(ns),length(phis),2);
CIc = CIr;
scl = norminv(0.975);

for in = 1:length(ns)
    n = ns(in);
    for iphi = 1:length(phis)
        phi = phis(iphi);
        
        %Generate panel data
        X = zeros(T,n);
        alpha = normrnd(0,1,1,n);
        ind = alpha>0.6;
        X(3:4,ind) = 1;
        epsilon = zeros(T+1,n);
        epsilon(1,:) = normrnd(0,1,1,n);
        u = normrnd(0,1,T,n);
        for tt = 1:T
            epsilon(tt+1,:) = phi*epsilon(tt,:) + u(tt,:);
        end
        Y = X*beta0 + repmat(pdelta,1,n) + repmat(alpha,T,1) + epsilon(2:end,:);
        t = repmat(diag([0 1 1 1]),n,1);
        t=t(:,2:end);
        
        %OLS
        b_ols(in,iphi,:) = ([ones(T*n,1) X(:) t])\Y(:);
        
        %FE
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([dX(:) dt]);
        b_fe(in,iphi,:) = XX\dY(:);
        
        %Standard errors
        r = dY(:) - XX*reshape(b_fe(in,iphi,:),4,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(4,4);
        for i = 1:n
            clustsum = clustsum + XX(T*(i-1)+1:T*i,:)' *r(T*(i-1)+1:T*i,:) * r(T*(i-1)+1:T*i,:)' * XX(T*(i-1)+1:T*i,:);
        end
        Vc = inv(XX'*XX)*clustsum*inv(XX'*XX);
        
        %Confidence intervals
        CIr(in,iphi,1) = b_fe(in,iphi,1) - scl*sqrt(Vr(1));
        CIr(in,iphi,2) = b_fe(in,iphi,1) + scl*sqrt(Vr(1));
        CIc(in,iphi,1) = b_fe(in,iphi,1) - scl*sqrt(Vc(1));
        CIc(in,iphi,2) = b_fe(in,iphi,1) + scl*sqrt(Vc(1));
    end
end

%% Simulation
nsim = 10000;
b_ols_sim = zeros(length(ns),length(phis),5,nsim);
b_fe_sim = zeros(length(ns),length(phis),4,nsim);
CRr = zeros(length(ns),length(phis),nsim);
CRc = CRr;

for isim = 1:nsim
for in = 1:length(ns)
    n = ns(in);
    for iphi = 1:length(phis)
        phi = phis(iphi);
        
        %Generate panel data
        X = zeros(T,n);
        alpha = normrnd(0,1,1,n);
        ind = alpha>0.6;
        X(3:4,ind) = 1;
        epsilon = zeros(T+1,n);
        epsilon(1,:) = normrnd(0,1,1,n);
        u = normrnd(0,1,T,n);
        for tt = 1:T
            epsilon(tt+1,:) = phi*epsilon(tt,:) + u(tt,:);
        end
        Y = X*beta0 + repmat(pdelta,1,n) + repmat(alpha,T,1) + epsilon(2:end,:);
        
        t = repmat(diag([0 1 1 1]),n,1);
        t = t(:,2:end);
        
        %OLS
        b_ols_sim(in,iphi,:,isim) = ([ones(T*n,1) X(:) t])\Y(:);     
        
        %FE
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([dX(:) dt]);
        b_fe_sim(in,iphi,:,isim)= XX\dY(:);
        
        %Standard errors
        r = dY(:) - XX*reshape(b_fe_sim(in,iphi,:,isim),4,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(4,4);
        for i = 1:n
            clustsum = clustsum + XX(T*(i-1)+1:T*i,:)' *r(T*(i-1)+1:T*i,:) * r(T*(i-1)+1:T*i,:)' * XX(T*(i-1)+1:T*i,:);
        end
        Vc = inv(XX'*XX)*clustsum*inv(XX'*XX);
        
        %Confidence intervals
        CIr_sim(1) = b_fe_sim(in,iphi,1,isim) - scl*sqrt(Vr(1));
        CIr_sim(2) = b_fe_sim(in,iphi,1,isim) + scl*sqrt(Vr(1));
        CIc_sim(1) = b_fe_sim(in,iphi,1,isim) - scl*sqrt(Vc(1));
        CIc_sim(2) = b_fe_sim(in,iphi,1,isim) + scl*sqrt(Vc(1));
        
        %Coverage rates
        CRr(in,iphi,isim) = (beta0>CIr_sim(1))*(beta0<CIr_sim(2));
        CRc(in,iphi,isim) = (beta0>CIc_sim(1))*(beta0<CIc_sim(2));    
    end
end
end

%Coverage Rates
CRr = mean(CRr,3);
CRc = mean(CRc,3);

%Average beta values
b_ols_sim = median(b_ols_sim,4);
b_fe_sim = median(b_fe_sim,4);

%% Exporting Results
mat = zeros(6,4);
mat2 = zeros(6,9);
for i=1:length(ns)
    for j=1:length(phis)
        row = (i-1)*2+j;
        mat(row,:) = [b_fe_sim(i,j,1) b_ols_sim(i,j,2) CRr(i,j) CRc(i,j)];
        mat2(row,:) = [b_fe(i,j,1) b_ols(i,j,2) reshape(b_fe(i,j,2:4),1,3) CIr(i,j,1) CIr(i,j,2) CIc(i,j,1) CIc(i,j,2)];
    end
end


tab = table(mat(:,1),mat(:,2),mat(:,3),mat(:,4), ...
    'RowNames',{'$\phi=0,n=40$' '$\phi=0.8,n=40$' '$\phi=0,n=70$' '$\phi=0.8,n=70$' '$\phi=0,n=100$' '$\phi=0.8,n=100$'}, ...
    'VariableNames',{'FE' 'OLS' 'Robust Coverage' 'Cluster Coverage'});
table2latex(tab,'ps6-1.tex')
tab2 = table(mat2(:,1),mat2(:,2),mat2(:,3),mat2(:,4),mat2(:,5),mat2(:,6),mat2(:,7),mat2(:,8),mat2(:,9), ...
    'RowNames',{'$\phi=0,n=40$' '$\phi=0.8,n=40$' '$\phi=0,n=70$' '$\phi=0.8,n=70$' '$\phi=0,n=100$' '$\phi=0.8,n=100$'}, ...
    'VariableNames',{'FE' 'OLS' '$\delta_2$' '$\delta_3$' '$\delta_4$' 'Robust LB' 'Robust UB' 'Cluster LB' 'Cluster UB'});
table2latex(tab2,'ps6-2.tex')