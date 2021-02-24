% pause('on')

clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Problem Sets\PS6'

%Written by Sarah Bass, 10/14/20

alpha = 0.01:0.01:0.49;

%Social Planner Equilibrium
s_l = (2*alpha+1)/3;
s_c = (1-alpha)/3;
s_g = (1-alpha)/3;
s_tau = s_g./(1-s_l);
s_u = log(s_l) + log(alpha + s_c) + log(alpha + s_g);

%Nash Equilibrium
n_l = alpha + 0.5;
n_c = 0.25 - 0.5*alpha;
n_g = 0.25 - 0.5*alpha;
n_tau = 0.5+0*n_g;
n_u = log(n_l) + log(alpha + n_c) + log(alpha + n_g);

%Ramsey Equilibrium
r_tau = 0*alpha;
objfn = @(t,a) -2*(1-t+a)^(-1) + 2*(1-t)^(-1) + (-3*a + 1 - 2*t)/(2*a*(1-t) + t*(1-t-a));
for i=1:length(alpha)
    tempobj = @(t)abs(objfn(t,alpha(i)));
    r_tau(i) = fmincon(tempobj,0.1,[1;-1],[1;0]);
end
r_l = 0.5+ (alpha./(2*(1-r_tau)));
r_c = (1-r_tau).*(1-r_l);
r_g = r_tau.*(1-r_l);
r_u = log(r_l) + log(alpha + r_c) + log(alpha + r_g);

%Figure for Ramsey Equilibrium
figure 
hold on;
plot (alpha, r_l, 'k')
plot (alpha, r_c, 'r')
plot (alpha, r_g, 'b')
plot (alpha, r_tau, 'g')
xlabel('\alpha')
legend ('Leisure', 'Consumption', 'Government Expenditure', 'Tax', 'Location', 'NorthEast')
title('Ramsey Equilibrium')
saveas(gcf,'ramsey.png')

%Figure comparing l across equilibriums
figure 
hold on;
plot (alpha, n_l, 'k')
plot (alpha, r_l, 'r')
plot (alpha, s_l, 'b')
xlabel('\alpha')
ylabel('Leisure')
legend ('Nash Equilibrium', 'Ramsey Equilibrium', 'Social Planners Equilibrium', 'Location', 'NorthEast')
title('Leisure Across Equilibria')
saveas(gcf,'leisure.png')

%Figure comparing c across equilibriums
figure 
hold on;
plot (alpha, n_c, 'k')
plot (alpha, r_c, 'r')
plot (alpha, s_c, 'b')
xlabel('\alpha')
ylabel('Consumption')
legend ('Nash Equilibrium', 'Ramsey Equilibrium', 'Social Planners Equilibrium', 'Location', 'NorthEast')
title('Consumption Across Equilibria')
saveas(gcf,'consumption.png')

%Figure comparing g across equilibriums
figure 
hold on;
plot (alpha, n_g, 'k')
plot (alpha, r_g, 'r')
plot (alpha, s_g, 'b')
xlabel('\alpha')
ylabel('Government Expenditure')
legend ('Nash Equilibrium', 'Ramsey Equilibrium', 'Social Planners Equilibrium', 'Location', 'NorthEast')
title('Government Expenditure Across Equilibria')
saveas(gcf,'government.png')

%Figure comparing u across equilibriums
figure 
hold on;
plot (alpha, n_u, 'k')
plot (alpha, r_u, 'r')
plot (alpha, s_u, 'b')
xlabel('\alpha')
ylabel('Utility')
legend ('Nash Equilibrium', 'Ramsey Equilibrium', 'Social Planners Equilibrium', 'Location', 'NorthEast')
title('Utility Across Equilibria')
saveas(gcf,'utility.png')

% Calculate beta
beta = 0*r_tau;
for i=1:length(alpha)
    if r_u(i)>n_u(i)
        temputil = log(r_l(i)) + log(alpha(i) + (0.5)*(1-r_l(i))) + log(alpha(i) + 0.5*(1-r_l(i)));
        genobj = @(b) abs(r_u(i)*(1/(1-b)) - temputil - (b/(1-b)*n_u(i)));
        beta(i) = fmincon(genobj,0.1,[1;-1],[1;0]);
    else
        beta(i) = 1;
    end
end

%Figure comparing u across equilibriums
figure
plot(alpha,beta, 'k')
title('\beta_L')
xlabel('\alpha')
saveas(gcf,'beta.png')