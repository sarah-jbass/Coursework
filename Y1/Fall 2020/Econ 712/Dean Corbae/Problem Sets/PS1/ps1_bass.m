% PROGRAM NAME: ps1_bass
%
% This program generates and plots the price dynamics given the first-order 
% difference equation discussed in the problem set given some initial price 
%
% Prepared by Fu Tan, modified by Sarah Bass
% Last updated: 09/09/15
pause('on')

clear;
clc;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Problem Sets\PS1'

%%Question 3

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
r = 0.01; % interest rate
d = 1; % constant dividend
p0a = 100; % initial price of 100
p0b = 90; % initial price of 90
p0c = 110; % initial price of 110
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
pvector_a = zeros(dim+1,1);
pvector_b = zeros(dim+1,1);
pvector_c = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99
pvector_a(1) = p0a; % giving value to the first element of the price vector
                % with the initial price
pvector_b(1) = p0b;
pvector_c(1) = p0c;

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%
for n = 2:dim+1 % starting from t = 1 to t = 100 
    pvector_a(n) = (1+r)*pvector_a(n-1)-d; % updating the price in the next period 
    pvector_b(n) = (1+r)*pvector_b(n-1)-d; % with the first-order difference
    pvector_c(n) = (1+r)*pvector_c(n-1)-d; % equation
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
hold on; 
plot(tvector(:),pvector_b(:));
plot(tvector(:),pvector_a(:));
plot(tvector(:),pvector_c(:));
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_{0a} = 90','P_{0b} = 100','P_{0c} = 110','Location','Northeast');
axis([0 dim 50 150])
saveas(gcf,'Question 3.png')

%%Question 2

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
r_pd = .1; % Slope is too close to 1 for the phase diagram to look good. 
        % I am using a different value for r so it looks better. 
x = 0:20; %just a grid for the phase diagram
xp = x*(1+r_pd)-d; %pt+1 = (1+r)pt - d
xss = 10+0*x; %steady state value

%%%%%%%%
% PLOTS
%%%%%%%%
figure
hold on
plot(x,x,'-.')
plot(x,xp)
plot(x,xss,'k--')
xline(xss(1),'k--')
hold off
legend('p_{t+1} = p_{t}','p_{t+1} =(1+r) p_t - d','p^{*}','Location','Northeast')
xlabel('p(t)')
xticks([10])
xticklabels({'p^*'})
ylabel('p(t+1)')
yticks([10])
yticklabels({'p^*'})
title('Phase Diagram')
saveas(gcf,'Question 2 Phase Diagram.png')

figure
hold on
plot(tvector(:),pvector_a)
plot(tvector(:),pvector_b)
plot(tvector(:),pvector_c)
hold off
legend('p_0 = p^{*}','p_0 < p^{*}','p_0 > p^{*}','Location','Northeast')
axis([0 dim 50 150])
xlabel('Time t'); ylabel('Price P_t');
title('P_t over time')
saveas(gcf,'Question 2 Pt over time.png')

%%Question 4

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
announ = 20;
new = 50; % point in time at which fed changes interest rates
rn = 0.02; % new interest rate

pvector_a(new) = 50; % set price vector at t=50, new p^* = 1/0.02

for n = 2:announ % from t = 1 to t = 20
    pvector_a(n) = (1+r)*pvector_a(n-1)-d; % updating the price in the next period 
                                       % with the first-order difference
                                       % equation
end

for n = new+1:dim+1 % from t = 51 to t = 100
    pvector_a(n) = (1+rn)*pvector_a(n-1)-d; % using new interest rate
end

for n = new-1:-1:announ+1 % from t = 49 to t = 21, backwards
   pvector_a(n)=(pvector_a(n+1)+d)/(1+r); % rearrange equation in terms of
                                          % n and n+1
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector_a(:));
title('Price Dynamics: Federal Funds Rate Increase');
xlabel('Time t'); ylabel('Price P_t');
legend('Price over time','Location','Northeast');
ylim([40 110])
saveas(gcf,'Question 4.png')
