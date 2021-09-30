clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Noah Williams\Problem Sets\PS2'

%Written by Sarah Bass, 11/18/20

%Parameters
beta = 0.95;
delta = 0.1;
z0 = 1;
z1 = 1.2;
gamma0 = 2;
gamma1 = 1.01;
psi = 0.35;

%Grid
k = 0.001:0.01:5; 

[kss0,css0,c,kprime,v] = iter_fn(k,beta,delta,z0,gamma0,psi);

%Phase lines
delk0 = z0 * k.^(psi) - delta*k; 
delc0 = z0 * k.^psi + (1-delta)*k - (((1/beta) -(1-delta))/(0.35*z0))^(1/(psi-1)) ;

%Figures
figure
hold on
plot(k,c,'k--')
plot(k,delk0,'b-')
plot(k,delc0,'r-')
plot(kss0,css0,'ko')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','Location','SouthEast')
title('3A Phase Diagram')
xlabel('k')
ylabel('c')
annotation('arrow',[0.2 0.2 ],[0.5 0.55])
annotation('arrow',[0.2 0.25 ],[0.5 0.5])
annotation('arrow',[0.87 0.87 ],[0.8 0.75])
annotation('arrow',[0.87 0.82 ],[0.8 0.8])
annotation('arrow',[0.5 0.5 ],[0.85 0.9])
annotation('arrow',[0.5 0.45 ],[0.85 0.85])
annotation('arrow',[0.65 0.65 ],[0.5 0.45])
annotation('arrow',[0.65 0.7 ],[0.5 0.5])
saveas(gcf,'3A.png')

figure
subplot(2,1,1)
plot(k,v)
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(2,1,2)
plot(k,kprime)
title('k-prime from value function')
ylabel('k-prime(k)')
set(gcf,'Color',[1 1 1])
saveas(gcf,'3A-2.png')

delk1 = z0 * k.^(psi) - delta*k; %phase
delc1 = z0 * k.^psi + (1-delta)*k - (((1/beta) -(1-delta))/(0.35*z0))^(1/(psi-1)) ;

[kss1,css1,c1,kprime1,v1] = iter_fn(k,beta,delta,z0,gamma1,psi);
figure
hold on
plot(k,c,'k')
plot(k,delk0,'r-')
plot(k,delc0,'b-')
plot(kss0,css0,'ko')
plot(k,c1,'k--')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','new saddle path (c(k))','Location','SouthEast')
title('3B Saddle Path')
xlabel('k')
ylabel('c')
saveas(gcf,'3B.png')

figure
subplot(3,1,1)
hold on
plot(k,v)
hold off
legend('old','Location','SouthEast')
title('Value function V(k)')
subplot(3,1,2)
plot(k,v1)
hold on
hold off
legend('new','Location','SouthEast')
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(3,1,3)
plot(k,kprime)
hold on
plot(k,kprime1)
hold off
legend('old','new','Location','SouthEast')
title('k-prime from value function')
ylabel('k-prime(k)')
set(gcf,'Color',[1 1 1])
saveas(gcf,'3B-2.png')

delk2 = z1 * k.^(psi) - delta*k; %phase
delc2 = z1 * k.^psi + (1-delta)*k - (((1/beta) -(1-delta))/(0.35*z1))^(1/(psi-1)) ;

[kss2,css2,c2,kprime2,v2] = iter_fn(k,beta,delta,z1,gamma0,psi);
figure
hold on
plot(k,c,'k')
plot(k,delk0,'r-')
plot(k,delc0,'b-')
plot(kss0,css0,'ko')
plot(k,c2,'k--')
plot(k,delk2,'r--')
plot(k,delc2,'b--')
plot(kss2,css2,'kx')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','new saddle path (c(k))','new \Delta k = 0','new \Delta c = 0','new SS','Location','SouthEast')
title('3C Saddle Path')
xlabel('k')
ylabel('c')
saveas(gcf,'3C.png')

figure
subplot(2,1,1)
hold on
plot(k,v)
plot(k,v2)
hold off
legend('old','new')
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(2,1,2)
plot(k,kprime)
hold on
plot(k,kprime2)
hold off
legend('old','new','Location','SouthEast')
title('k-prime from value function')
ylabel('k-prime(k)')
set(gcf,'Color',[1 1 1])
saveas(gcf,'3C-2.png')

nt=50;
%find saddle path for transition dynamics
cc0 = css0:0.00001:css2;
kk = 0*cc0+kss0;
cc=cc0;
for i=1:nt
    [kk,cc] = LOM(kk,cc,z1,psi,delta,beta,gamma0);
end
[~,ind] = min(abs(cc - css2)+abs(kk-kss2));
cc0 = cc0(ind);
ctraj = cc0*ones(1,nt);
ktraj = kss0*ones(1,nt);
for i=2:nt
    [ktraj(i),ctraj(i)] = LOM(ktraj(i-1),ctraj(i-1),z1,psi,delta,beta,gamma0);
end

ct = [css0 ctraj css2];
kt = [kss0 ktraj kss2];
figure
hold on
plot(kt,ct,'k--')
plot(kt(1),ct(1),'ko')
plot(kt(2),ct(2),'k+')
plot(kt(end),ct(end),'kx')
hold off
set(gcf,'Color',[1 1 1])
title('Computational solution to particle trajectory for part (c)')
xlabel('k')
ylabel('c')
legend('Particle trajectory','initial steady state','jump','final steady state','Location','SouthEast')
saveas(gcf,'traj.png')

%Functions
function [kk,cc] = LOM(k,c,z,psi,delta,beta,gamma)
    kk = z*k.^psi - c + (1-delta)*k;
    cc = c.*(beta*(psi*z*k.^(psi-1)+1-delta)).^(1/gamma);
end
    
function [kss,css,c,kprime,v] = iter_fn(k,beta,delta,z0,gamma0,psi)
    % Iterative method of solving value function

    kss = (((1/beta) - 1 + delta)/(psi*z0))^(1/(psi - 1));
    css = z0*kss^(psi) - delta * kss;

    nk=length(k);
    v0 = (1/((1-gamma0)*(1-beta)))*((1-beta)*k).^(1-gamma0); % initial guess of value function approx geometric
    diff = 1e6;
    maxiter = 1e6;
    tol=1e-6;
    iter=1;
    posval = ones(nk); %prealloc
    
    for i=1:nk % defines what is possible
    posval(i,k>z0*k(i).^(psi) + (1-delta)*k(i)) = 0; %k' legal: c=0: max is pz*k^(ppsi) + (1-pdelta)*k
    end
    
    posval = logical(posval);
    v=0*v0;
    kprime = 0*v0;
    
    % iterate until converged
    while (diff>tol&&iter<maxiter)
        if mod(iter,100)==0; disp(['iteration: ' num2str(iter) ', diff = ' num2str(diff) ]); end
        for i=1:nk
            %k0=k(i); %If I start the period with k0 capital
            kpopt = k(posval(i,:)); % I have these options for k'
            vi = (z0*k(i).^(psi) + (1-delta)*k(i) - kpopt).^(1-gamma0)./(1-gamma0) + beta*v0(posval(i,:));   %which yield me this much utility
            [v(i),ind] = max(vi); % value function is the best option
            kprime(i) = kpopt(ind);
        end
        diff = norm(v-v0);
        iter=iter+1;
        v0=v;
    end
    c = z0*k.^(psi) + (1-delta)*k - kprime;
end
