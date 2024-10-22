clear
clc

% Sample Solution code of Problem Set 2, ECON712
% Written by Anson Zhou

%% Problem 1
% Parameter setup
z = 1;
alpha = 0.3;
delta = 0.1;
beta = 0.97;

% steady state
k_bar = % fill in the blank for steady state values; 
c_bar =  ;

% compute the Jacobian matrix at steady state
J_11 =  ;  % fill in the blank for elements of J
J_12 =  ;
J_21 =  ;
J_22 =  ;
J = [J_11, J_12; J_21, J_22];

% eigenvalues and eigenvectors
[V,D]=eig(J);

% new steady state
z = z + 0.1;
k_bar2 = % fill in the blank, calculate new steady state values;
c_bar2 = ;

% Jacobian matrix at new steady state
J2_11 = % new Jacobi matrix ;
J2_12 = ;
J2_21 = ;
J2_22 = ;
J2 = [J2_11, J2_12; J2_21, J2_22];

% eigenvalues and eigenvectors
[V2,D2]=eig(J2);


%% determine coefficients for specific solution
t0 = 5;
m2 = % calculate the constant needed for the particular solution, I call it m2;

% consumption response at t0
c_t0 = c_bar2 + V2(2,2)*m2*(D2(2,2).^t0);

% trajectory using linearized saddle path
Prd = 20-t0+1;  % number of periods we need to calculate

% Hereafter I use Trj (a 2 by Prd matrix) to store the values of capital
% and consumption for the transition period. First row is capital, second
% row is consumption
Trj1 = zeros(2,Prd);
Trj1(1,1) = k_bar;   % At t0, agent are stuck with k_bar
Trj1(2,1) = c_t0;     % They can choose consumption so that they will fall onto the new saddle path
for i = 2:Prd
    Trj1(:,i)= % Fill in the system dynamics, i.e. relating Trj1(:,i) to Trj1(:,i-1) using 2(c);
end


% trajectory using original equation and the "jump" of c_5 calculated under
% linearization 
Trj2 = zeros(2,Prd);
Trj2(1,1)=k_bar;
Trj2(2,1)=c_t0;
for i = 2:Prd
    k = Trj2(1,i-1);
    c = Trj2(2,i-1);
    Trj2(1,i) = % Fill out this line using the actual nonlinear dynamics (Hint: equation 1(a));
    Trj2(2,i) = % Hint: 1(b);
end

% The matrix "before" is simply storing the system values up to date t0, it
% is used to graph the whole system from date 0 to date 20
before = repmat([k_bar; c_bar],1, t0-1);
Trj1_whole = [before, Trj1];   % Transition under linear hypothesis
Trj2_whole = [before, Trj2];   % Actual transition using "jump" from linearization

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

% This part finds the actual c_t0 needed using "shooting method". The idea
% is to try different values of c_t0 and choose the one such that the
% system will converge to the new steady state after a long period of time

M = 10000;  % Picking M number of candidates from c_t0
C_range = linspace(c_bar, c_bar2, M);  % Choosing candidates of c_t0
Distance = zeros(1,M);  % The matrix that stores the distance between the new steady state
% and where the system is at after a long period of time following
% particular c_t0
for i = 1:M
    Distance(1,i)=calib(% Read "calib.m" and fill in the blanks, you might need to change the file name);
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
    Trj3(1,i) = % Fill in the blank;
    Trj3(2,i) = % Fill in the blank;
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

