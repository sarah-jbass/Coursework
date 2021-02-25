clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 714\Dmitry Mukhin\Problem Sets\PS4'

%Written by Sarah Bass, 2/20/21

%Params
K = 100000;
Nk = 20;
Aik= exp(normrnd(0,1,K,Nk));
W=1;
epsilon = 1e-5;
prho = 1+ epsilon;
ptheta = 5;

s0 = (1/Nk)+ 0*Aik;
s = s0;
diff = 999;
tol = 1e-5;
maxiter = 1e10;
iter = 1;
tune = 0.6;

while (iter<maxiter) && (diff>tol)
    pik = (W./Aik).*(1-(1./((1-ptheta) + s.*(ptheta - prho))));
    pk = (sum(pik.^(1-ptheta),2)).^(1/(1-ptheta));
    s0 = s;
    s = (pik./pk).^(1-ptheta);
    diff = sum(sum(abs(s-s0),2));
    s = s0*tune + s*(1-tune);
    iter = iter + 1;
end

p = ((1/K)*sum(pk.^(1-prho))).^(1/(1-prho));
c = W/p;


% Initial guess of s
% Given s, calculate p_ik
% Given p_ik, calculate p_k
% Given p_k, P_ik, calculate s
% Check how close s is to s0
% New s is a weighted avg using a tuning parameter (0.6)
% Calculate p, c = w/p
