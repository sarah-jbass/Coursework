clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Noah Williams\Problem Sets\PS3'

%Written by Sarah Bass, 12/1/20

%Part B
Q = [0.85 0.15; 0.05 0.95]; %Q
P0 = [1 0]; %P
tol = 1e-5; %tolerance
maxiter = 1e5; %maximum number of iterations
iter=1; %first iteration
diff = 1; %starting difference, doesn't matter
while ((iter<maxiter)&&(diff>tol))
    P=P0*Q;
    diff = norm(P-P0); %use norm bc vector, finding distance between p and p0
    P0 = P;
    iter=iter+1;
end

%Part C
pbeta = 0.95; %beta
pgamma = 3; %gamma
pr = 0.03; %interest rate
pw = 1.1; %wage
l = [0.7 1.1]; %labor endowment
assets = 0:1e-2:3; % asset grid
na = length(assets); %number of assets
V0 = [log(assets') log(assets')]; %initial guess for value function
legal = true(na,na,2); %identify assets that are feasible
for nl = 1:2
   for nap = 1:na
       legal(nap,:,nl) = (pw*l(nl) + (1+pr)*assets(nap) - assets')>0;
   end
end

V=0*V0; %preallocation, value doesn't matter but same size as V
aprime = V; %preallocation
aind = V; %preallocation
diff=999; %reset difference
iter=1; %reset iterations
while ((iter<maxiter)&&(diff>tol))
    disp(['beginning iteration ' num2str(iter) ', diff = ' num2str(diff)])
    for nl = 1:2 %each labor endowment
    for ai = 1:na %and each asset holding
        Val = (pw*l(nl) + (1+pr)*assets(ai) - assets(legal(ai,:,nl))' ).^(1-pgamma)./(1-pgamma) ...
            + pbeta*(V0(legal(ai,:,nl),1)*Q(nl,1) + V0(legal(ai,:,nl),2)*Q(nl,2));
        [V(ai,nl),ind] = max(Val); %fill in V
        aprime(ai,nl) = assets(ind); %fill in assets w/ index for max(Val)
        aind(ai,nl) = ind; %fill in index
    end
    end
    iter = iter+1;
    diff = sum(abs(V(:)-V0(:))); %slower rate of convergence than norm, functionally the same
    V0 = V;
end
% save results
% else
%     load results
% end

abar = zeros(2,1);
for nl = 1:2
    %for ai= 1:na
        [~,ind]=min(abs(assets' - aprime(:,nl))); %takes second argument of min function, which is the index
        abar(nl) = ind(1);
        % x: assets(abar(1))
        % y: aprime(abar(1))
    %end
end

%Value function plot
figure
hold on
plot(assets,V(:,1),'k')
plot(assets,V(:,2),'r')
hold off
title('Value function')
xlabel('Assets')
ylabel('Value')
legend('Low labor','High labor','Location','SouthEast')
set(gcf,'Color',[1 1 1])
xlim([-0.1 3])
saveas(gcf,'value.png')

%Next period assets plot
figure
hold on
plot(assets,aprime(:,1),'k')
plot(assets,aprime(:,2),'r') %plot(assets,aprime(aprime(:,1)<assets,1))
plot(assets,assets,'b--')
plot(assets(abar(1)), aprime(abar(1),1), 'rx')
plot(assets(abar(2)), aprime(abar(2),2), 'rx')
hold off
title('Asset holdings for next period')
xlabel('Assets')
ylabel('a prime')
legend('Low labor','High labor', 'a-prime = a', 'Cutoff', 'Location','SouthEast')
set(gcf,'Color',[1 1 1])
xlim([-0.1 3])
saveas(gcf,'assets.png')

%Part D

% Stationary distribution
dist0 = ones(size(V));
dist0 = dist0./(sum(dist0(:))); %every weight is 1/2na
dist = dist0;
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    dist = 0*dist0;
    for nl = 1:2
       for ai = 1:na
           if dist0(ai,nl)>1e-14 % 
               targ = aind(ai,nl);%assets==aprime(ai,nl);
               dist(targ,1) = dist(targ,1)+dist0(ai,nl)*Q(nl,1);
               dist(targ,2) = dist(targ,2)+dist0(ai,nl)*Q(nl,2);
           end
       end
    end
    iter = iter+1;
    diff = sum(abs((dist0(:)-dist(:))));
    dist0 = dist;
end

%Checks!
sum(dist(:,1)) %should =.25
sum(dist(:,2)) %should = .75
sum(dist(:)) %should =1

%Marginal distribution graph
figure
hold on
plot(assets,dist(:,1),'k')
plot(assets,dist(:,2),'r')
hold off
title('Marginal asset distribution')
xlabel('Assets')
ylabel('Probability mass')
legend('Low Labor', 'High Labor', 'Location','NorthEast')
set(gcf,'Color',[1 1 1])
xlim([-0.1 3])
saveas(gcf,'marginal.png')