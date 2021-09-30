% pause('on')

clear;
clc;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Problem Sets\PS2'

% Sample Solution code of Problem Set 2, ECON712
% Written by Anson Zhou, edited by Sarah Bass

%% Problem 1
% Parameter setup
pz = 1;
palpha = 0.3;
pdelta = 0.1;
pbeta = 0.97;

% steady state
k_bar = ((pz*pbeta*palpha)/(1-(pbeta*(1-pdelta))))^(-1/(palpha-1)) ; 
c_bar =  pz*k_bar^palpha - pdelta*k_bar;

% compute the Jacobian matrix at steady state
syms k c z alpha delta beta
J = jacobian([z*k^alpha + (1-delta)*k - c, c*beta*(1-delta+alpha*z*(z*k^alpha + (1-delta)*k - c)^(alpha-1))],[k, c]);

varnames = {'k', 'c', 'z', 'alpha', 'delta', 'beta'};
params = {'k_bar', 'c_bar', 'pz', 'palpha', 'pdelta', 'pbeta'};

J2 = J %Creating a Jacobian matrix w/ the param values inserted, solved at steady state
for i = 1:length(varnames) 
   eval(['J2=subs(J2,' varnames{i} ',' params{i} ');']);
end 
J2 = double(J2)

[V,D] = eig(J2) %These are the eigenvalues and eigenvectors for the original Jacobian at steady state
slope = V(2,2)/V(1,2)

% new steady state
pz = pz + 0.1;
k_bar2 = ((pz*pbeta*palpha)/(1-(pbeta*(1-pdelta))))^(-1/(palpha-1));
c_bar2 = pz*k_bar2^palpha - pdelta*k_bar2;

params2 = {'k_bar2', 'c_bar2', 'pz', 'palpha', 'pdelta', 'pbeta'};

J3 = J %Creating a Jacobian matrix w/ the new param values inserted, solved at new steady state
for i = 1:length(varnames) 
   eval(['J3=subs(J3,' varnames{i} ',' params2{i} ');']);
end 
J3 = double(J3)

[V2,D2] = eig(J3) %These are the eigenvalues and eigenvectors for the new Jacobian at steady state

% %% determine coefficients for specific solution
t0 = 5;
c2 = (k_bar - k_bar2)/(V2(1,2)*(D2(2,2).^t0)) % calculate the constant needed for the particular solution, I call it m2;

% consumption response at t0
c_t0 = c_bar2 + V2(2,2)*c2*(D2(2,2).^t0);

% trajectory using linearized saddle path
Prd = 20-t0+1;  % number of periods we need to calculate

% Hereafter I use Trj (a 2 by Prd matrix) to store the values of capital
% and consumption for the transition period. First row is capital, second
% row is consumption
Trj1 = zeros(2,Prd);
Trj1(1,1) = k_bar;   % At t0, agent are stuck with k_bar
Trj1(2,1) = c_t0;     % They can choose consumption so that they will fall onto the new saddle path
for i = 2:Prd
    Trj1(1,i)= V2(1,2)*c2*(D2(2,2).^(t0+i-1))+k_bar2;
    Trj1(2,i)= V2(2,2)*c2*(D2(2,2).^(t0+i-1))+c_bar2; % Fill in the system dynamics, i.e. relating Trj1(:,i) to Trj1(:,i-1) using 2(c);
end

% trajectory using original equation and the "jump" of c_5 calculated under
% linearization 
Trj2 = zeros(2,Prd);
Trj2(1,1)=k_bar;
Trj2(2,1)=c_t0;
for i = 2:Prd
    k = Trj2(1,i-1);
    c = Trj2(2,i-1);
    Trj2(1,i) = pz*k^palpha + (1-pdelta)*k - c; % Fill out this line using the actual nonlinear dynamics (Hint: equation 1(a));
    Trj2(2,i) = pbeta*c*(1 - pdelta + palpha*pz*(pz*k^palpha + (1 - pdelta)*k - c)^(palpha - 1)); % Hint: 1(b);
end
% 
% % The matrix "before" is simply storing the system values up to date t0, it
% % is used to graph the whole system from date 0 to date 20
before = repmat([k_bar; c_bar],1, t0-1);
Trj1_whole = [before, Trj1];   % Transition under linear hypothesis
Trj2_whole = [before, Trj2];   % Actual transition using "jump" from linearization
% 
time = 1:20;
 
% Linear plot
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(2,:),'-*')
xlabel('Time')
title('c_t')
saveas(gcf,'Question1_1.png')

% Showing linear value under actual dynamics does not converge
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(2,:),'-*')
xlabel('Time')
title('c_t')
saveas(gcf,'Question1_2.png')

% This part finds the actual c_t0 needed using "shooting method". The idea
% is to try different values of c_t0 and choose the one such that the
% system will converge to the new steady state after a long period of time

M = 10000;  % Picking M number of candidates from c_t0
C_range = linspace(c_bar, c_bar2, M);  % Choosing candidates of c_t0
Distance = zeros(1,M);  % The matrix that stores the distance between the new steady state
% and where the system is at after a long period of time following
% particular c_t0
for i = 1:M
    c_in = C_range(1,i);
    Distance(1,i)=calib_Bass(c_in, palpha, pbeta, pdelta, pz, k_bar, k_bar2, c_bar2);% Read "calib.m" and fill in the blanks, you might need to change the file name);
    % The calib.m provided is a function that calibrates where the system
    % is at after 50 periods
end
[DD,I] = min(Distance);     % Picking the particular c_t0 that minimizes the distance
c_star = C_range(I);

% Calibrating the actual nonlinear transition using c_star chosen above.
Trj3 = zeros(2,Prd);
Trj3(1,1)=k_bar;
Trj3(2,1)=c_star;
for i = 2:Prd
    k = Trj3(1,i-1);
    c = Trj3(2,i-1);
    Trj3(1,i) = pz*k^palpha + (1-pdelta)*k - c;
    Trj3(2,i) = pbeta*c*(1 - pdelta + palpha*pz*(pz*k^palpha + (1 - pdelta)*k - c)^(palpha - 1));
end

Trj3_whole = [before, Trj3];  

% Actual transition using "shooting algorithm"
% Actual transition, using c_5 calculated from "shooting algorithm"
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(2,:),'-*')
xlabel('Time')
title('c_t')
saveas(gcf,'Question1_3.png')
