clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Noah Williams\Problem Sets\PS4'

%Written by Sarah Bass, 12/7/20

%Params
pbeta = 0.95;
pdelta = 0.1;
pz0 = 1;
pgamma0 = 2;
pgamma1 = 1.01;
palpha = 0.35;

%% 2A - optimal policy function
rerun = 0;
if rerun == 1
    k = 0.001:0.001:5; %k grid
    [c0,kprime0,v0] = iter(k,pbeta,pdelta,pz0,pgamma0,palpha);
    z = [0.8 1.2];
    P = [0.9 0.1;0.1 0.9];
    [c,kprime,v] = iter2(k,pbeta,pdelta,z,P,pgamma0,palpha);
    save res_q2
else
    load res_q2
end

%Plots
figure
hold on
plot(k,kprime(:,1),'r')
plot(k,kprime(:,2),'b')
plot(k,kprime0,'k')
hold off
set(gcf,'Color',[1 1 1])
legend('low productivity','high productivity','deterministic','Location','SouthEast')
title('k prime vs k')
ylabel('k prime')
xlabel('k')
saveas(gcf,'q2capital.png')

figure
hold on
plot(k,c(:,1),'r')
plot(k,c(:,2),'b')
plot(k,c0,'k')
hold off
set(gcf,'Color',[1 1 1])
legend('low productivity','high productivity','deterministic','Location','SouthEast')
title('c vs k')
ylabel('c')
xlabel('k')
saveas(gcf,'q2consumption.png')


%% 2B - simulate
rerun2=0;
if rerun2==1
T = 1e5;
burn = floor(T/10);
dist0 = ones(size(k'));
dist0 = dist0./sum(dist0(:));
stat = 1; % initial state
cons = zeros(T,1);
Y = cons;
kap = cons;
I = cons;
w = cons;
r = cons;
% make transition matrices
mat1 = zeros(length(k));
mat2 = mat1;
for i=1:length(k)
    mat1(:,i) = 0+kprime(i,1)==k;
    mat2(:,i) = 0+kprime(i,2)==k;
end

for tt = 1:T+burn
    if mod(tt,burn)==0; disp(['iteration ' num2str(tt)]); end
    draw = rand;
    if draw<0.1; stat = 1-stat; end % z transition with probability 0.1
    switch stat
        case 0
            dist = mat1*dist0;
        otherwise
            dist = mat2*dist0;
    end
    if tt>burn % after burn-in, calculate 
        t = tt-burn;
        % calculate values 
        cons(t) = sum(c(:,stat+1).*dist0);
        kap(t) = sum(k'.*dist0);
        I(t) =sum(kprime(:,stat+1).*dist0 - (1-pdelta)*k'.*dist0);
        Y(t) = kap(t)^(palpha);
        r(t) = sum((palpha)*kap(t).^(palpha-1));
        w(t) = sum(kap(t).^(palpha) - r(t)*kap(t));
    end
    dist0=dist;
end
save simres
else
    load simres
end
% calculate deterministic steady state
matss = 0*mat1;
for i=1:length(k)
    matss(:,i) = 0+kprime0(i)==k;
end
init=ones(size(k'));
init=init./(sum(init));
tol=1e-6;
diff = 999;
maxiter=1e6;
iter=1;
while (iter<maxiter)&&(diff>tol)
    next = matss*init;
    iter = iter+1;
    diff = sum(abs(next - init));
    init = next;
end
consss = sum(c0(:).*next);
kapss = sum(k'.*next);
Yss = kapss^(palpha);
Iss =sum(kprime0'.*next - (1-pdelta)*k'.*next);
rss = sum((palpha)*kapss.^(palpha-1));
wss = sum(kapss.^(palpha) - rss*kapss);

consa = mean(cons);
kapa = mean(kap);
Ia = mean(I);
ra = mean(r);
wa = mean(w);
Ya = mean(Y);
% 
% Deterministic = [ kapss Yss consss Iss wss rss]';
% Simulation = [kapa Ya consa  Ia wa ra]';
% tab = table(Deterministic,Simulation,'RowNames',{'Capital' 'Output' 'Consumption' 'Investment' 'wages' 'interest rates'});
% table2latex(tab,'q2.tex')


%% 2C - volatility and correlation
cvol = std(cons);
Yvol = std(Y);
Ivol = std(I);
cycorr = corr(cons,Y);
rycorr = corr(r,Y);

% Simulation = [cvol Yvol Ivol cycorr rycorr]';
% tab = table(Simulation,'RowNames',{'Consumption Volatility' 'Output Volatility' 'Investment Volatility' 'C-Y Correlation' 'r-Y Correlation'});
% table2latex(tab,'q2vols.tex')
