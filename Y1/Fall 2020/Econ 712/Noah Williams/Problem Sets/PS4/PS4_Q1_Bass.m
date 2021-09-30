clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Noah Williams\Problem Sets\PS4'

%Written by Sarah Bass, 12/7/20

%Params
Q = [0.85 0.15; 0.05 0.95]; %Q
P0 = [1 0]; %P
tol = 1e-5; %tolerance
maxiter = 1e5; %maximum number of iterations
% iter=1; %first iteration
% diff = 1; %starting difference, doesn't matter

%Part C
palpha = 0.36; %alpha
pbeta = 0.95; %beta
pgamma = 3; %gamma,
pdelta = 0.08; %delta
% pr = 0.03; %interest rate
% pw = 1.1; %wage
l = [0.7 1.1]; %labor endowment
assets = 0:0.05:50; % asset grid
na = length(assets); %number of assets
V0 = [log(assets') log(assets')]; %initial guess for value function
kd0=5.0159;
tune = .95;

tol2 = 0.02;
maxiter2 = 1e5;

rerun = 0;
if rerun == 1
[kd] = find_kd(kd0,assets,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V0,tol2,maxiter2,tune);
[ks,V,c,aprime,dist] = find_ks(kd,assets,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V0);
save results
else
    load results
end

%Consumption CDF
cvec = sort(unique(c));
cpmf = 0*cvec;
for nl = 1:2
    for i = 1:na
        cpmf(cvec==c(i,nl)) = cpmf(cvec==c(i,nl))+dist(i,nl);
    end
end
ccdf = cumsum(cpmf);

%Income CDF
pr = (palpha)*(kd)^(palpha - 1); % FOC of production function
pw = kd^palpha - pr*kd; %zero profit condition
inc = pw*repmat(l,na,1) + pr*repmat(assets',1,2);
incvec = sort(unique(inc));
incpmf = 0*incvec;
for nl = 1:2
    for i= 1:na
        incpmf(incvec==inc(i,nl)) = incpmf(incvec==inc(i,nl))+dist(i,nl);
    end
end
inccdf = cumsum(incpmf);

%Consumption and income distribution plot
figure
subplot(2,1,1)
plot(cvec,ccdf,'k')
set(gcf,'Color',[1 1 1])
title('Cumulative Consumption Distribution')
xlabel('Consumption')
ylabel('CDF')
xlim([0 4])
subplot(2,1,2)
plot(incvec,inccdf,'k')
set(gcf,'Color',[1 1 1])
title('Cumulative Income Distribution')
xlabel('Income')
ylabel('CDF')
xlim([0 4])
saveas(gcf,'firstdist.png')

kd02 = 4.9721;
assets2 = -2:0.05:50; % asset grid
V02 = [log(assets2'+3) log(assets2'+3)]; %initial guess for value function
tune = 0.98;
rerun2=0;
if rerun2==1
[kd2] = find_kd(kd02,assets2,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V02,tol2,maxiter2,tune);
[ks2,V2,c2,aprime2,dist2] = find_ks(kd2,assets2,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V02);
save results2
else
    load results2
end

%New Consumption CDF
c2vec = sort(unique(c2));
c2pmf = 0*c2vec;
for nl = 1:2
    for i = 1:na
        c2pmf(c2vec==c2(i,nl)) = c2pmf(c2vec==c2(i,nl))+dist2(i,nl);
    end
end
c2cdf = cumsum(c2pmf);

%New Income CDF
na2 = length(assets2); %number of assets
pr2 = (palpha)*(kd2)^(palpha - 1); % FOC of production function
pw2 = kd2^palpha - pr2*kd2; %zero profit condition
inc2 = pw2*repmat(l,na2,1) + pr2*repmat(assets2',1,2);
inc2vec = sort(unique(inc2));
inc2pmf = 0*inc2vec;
for nl = 1:2
    for i= 1:na2
        inc2pmf(inc2vec==inc2(i,nl)) = inc2pmf(inc2vec==inc2(i,nl))+dist2(i,nl);
    end
end
inc2cdf = cumsum(inc2pmf);

%New consumption and income distribution plot
figure
subplot(2,1,1)
hold on
plot(cvec,ccdf,'k')
plot(c2vec,c2cdf,'r')
hold off
set(gcf,'Color',[1 1 1])
title('Cumulative Consumption Distribution')
xlabel('Consumption')
ylabel('CDF')
legend('Min asset = 0','Min asset = -2', 'Location','NorthEast')
xlim([0 4])
subplot(2,1,2)
hold on
plot(incvec,inccdf,'k')
plot(inc2vec,inc2cdf,'r')
hold off
set(gcf,'Color',[1 1 1])
title('Cumulative Income Distribution')
xlabel('Income')
ylabel('CDF')
legend('Min asset = 0','Min asset = -2', 'Location','NorthEast')
xlim([0 4])
saveas(gcf,'seconddist.png')
