% PROGRAM NAME: ps1
%
% This program generates and plots the price dynamics given the first-order 
% difference equation discussed in the problem set given some initial price 
%
% Prepared by Fu Tan
% Last update: 09/08/15

clear;
clc;

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
plot(tvector(:),pvector_a(:));
plot(tvector(:),pvector_b(:));
plot(tvector(:),pvector_c(:));
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_0a = 100','P_0b = 90','P_0c = 110','Location','Northeast');
axis([0 dim 50 150])
