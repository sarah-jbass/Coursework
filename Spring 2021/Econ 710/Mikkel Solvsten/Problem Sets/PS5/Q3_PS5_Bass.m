clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 710\Mikkel Solvsten\Problem Sets\PS5'

%Written by Sarah Bass, 2/23/21

%Params
rng(100); %seed
alpha0 = 100;
delta0 = 100;
beta0 = 1;
T = [50 150 250];
rho1 = [0.7 0.9 0.95];
betas = zeros(4,3,3);
CIl = zeros(4,3,3);
CIu = zeros(4,3,3);
scale = norminv(0.975); %width of CI, approx 1.96

%Question 2A
for i=1:3
    for j=1:3
        V = normrnd(0,1,T(i),1); 
        U = normrnd(0,1,T(i),1);
        Y = zeros(T(i)+1,1);
        X = zeros(T(i)+1,1);
        Y(1) = normrnd(0,1); 
        X(1) = normrnd(0,1); 
        for k=1:T(i)
            X(k+1) = 0.3*X(k) + V(k);
            Y(k+1) = alpha0 + k*beta0 + X(k+1)*delta0 + Y(k)*rho1(j) + U(k);
        end
        tX = [ones(T(i),1) (1:T(i))' X(2:end) Y(1:end-1)]; %regressors
        betas(:,i,j) = tX\Y(2:end); %solve system of eqns
        r = Y(2:end) - tX*betas(:,i,j);
        Vhc0 = inv(tX' * tX) * (tX' * diag(r.^2) * tX) * inv(tX' * tX);
        SE = sqrt(diag(Vhc0));
        CIl(:,i,j) = betas(:,i,j) - scale*SE;
        CIu(:,i,j) = betas(:,i,j) + scale*SE;
    end
end

%Question 2B
nsim = 10000;
betasN = zeros(4,3,3,nsim);
coverN = zeros(4,3,3,nsim);
for n = 1:nsim
    for i = 1:3
        for j = 1:3
            V = normrnd(0,1,T(i),1); 
        U = normrnd(0,1,T(i),1);
        Y = zeros(T(i)+1,1);
        X = zeros(T(i)+1,1);
        Y(1) = normrnd(0,1); 
        X(1) = normrnd(0,1); 
        for k=1:T(i)
            X(k+1) = 0.3*X(k) + V(k);
            Y(k+1) = alpha0 + k*beta0 + X(k+1)*delta0 + Y(k)*rho1(j) + U(k);
        end
        tX = [ones(T(i),1) (1:T(i))' X(2:end) Y(1:end-1)]; %regressors
        betasN(:,i,j,n) = tX\Y(2:end); %solve system of eqns
        r = Y(2:end) - tX*betasN(:,i,j,n);
        Vhc0 = inv(tX' * tX) * (tX' * diag(r.^2) * tX) * inv(tX' * tX);
        SE = sqrt(diag(Vhc0));
        Cl = betasN(:,i,j,n) - scale*SE;
        Cu = betasN(:,i,j,n) + scale*SE;
        if alpha0>Cl(1) && alpha0<Cu(1)
            coverN(1,i,j,n) = 1;
        end
        if beta0>Cl(2) && beta0<Cu(2)
            coverN(2,i,j,n) = 1;
        end
        if delta0>Cl(3) && delta0<Cu(3)
            coverN(3,i,j,n) = 1;
        end
        if rho1(j)>Cl(4) && rho1(j)<Cu(4)
            coverN(4,i,j,n) = 1;
        end
        end
    end
end

PcoverN = mean(coverN,4);
MbetasN = mean(betasN,4);

mats = zeros(9,9);
for i = 1:3
    for j=1:3
        mats(3*(i-1)+j,[1 4 7]) = betas([1 3 4],i,j)';
        mats(3*(i-1)+j,1+[1 4 7]) = CIl([1 3 4],i,j)';
        mats(3*(i-1)+j,2+[1 4 7]) = CIu([1 3 4],i,j)';
    end
end
tabs = table(mats(:,1),mats(:,2),mats(:,3),mats(:,4),mats(:,5),mats(:,6),mats(:,7),mats(:,8),mats(:,9),...
    'RowNames',{'$(T=50,\rho_1=0.7)$' '$(T=50,\rho_1=0.9)$' '$(T=50,\rho_1=0.95)$' ...
    '$(T=150,\rho_1=0.7)$' '$(T=150,\rho_1=0.9)$' '$(T=150,\rho_1=0.95)$' ...
    '$(T=250,\rho_1=0.7)$' '$(T=250,\rho_1=0.9)$' '$(T=250,\rho_1=0.95)$'}, ...
    'VariableNames',{'$\hat{\alpha}_0$' '$\hat{\alpha}_0$ LB' '$\hat{\alpha}_0$ UB' ...
    '$\hat{\delta}_0$' '$\hat{\delta}_0$ LB' '$\hat{\delta}_0$ UB'  ...
    '$\hat{\rho}_1$' '$\hat{\rho}_1$ LB' '$\hat{\rho}_1$ UB'});
table2latex(tabs,'ps5s.tex')
% construct table for simulation results
mat = zeros(9,6);
for i = 1:3
    for j=1:3
        mat(3*(i-1)+j,[1 3 5]) = MbetasN([1 3 4],i,j)';
        mat(3*(i-1)+j,[2 4 6]) = PcoverN([1 3 4],i,j)'; 
    end
end
tab = table(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),...
    'RowNames',{'$(T=50,\rho_1=0.7)$' '$(T=50,\rho_1=0.9)$' '$(T=50,\rho_1=0.95)$' ...
    '$(T=150,\rho_1=0.7)$' '$(T=150,\rho_1=0.9)$' '$(T=150,\rho_1=0.95)$' ...
    '$(T=250,\rho_1=0.7)$' '$(T=250,\rho_1=0.9)$' '$(T=250,\rho_1=0.95)$'}, ...
    'VariableNames',{'$\hat{\alpha}_0$ Mean' '$\hat{\alpha}_0$ Coverage' ...
    '$\hat{\delta}_0$ Mean' '$\hat{\delta}_0$ Coverage' '$\hat{\rho}_1$ Mean' '$\hat{\rho}_1$ Coverage'});
table2latex(tab,'ps5.tex')
        
        