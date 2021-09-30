

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section estimates the Berry model, but ignores the order of entry effects.
%Although we obtain a unique estimate of the number of firms, N, in the market, the 
%estimates do not pertain to a unique equilibrium.  It is shown in Berry's paper that
%an order of entry is imposed then the estimates pertain to a unique equilibrium solution.



%Load the data and specify market data x, firm data xi, instruments z, and 
%market entry information y.  I believe all the instruments are the same, except I
%excluded
clear all
run loadWideData;
const=ones(size(pop,1),1);
x=[const pop distance];
xi=[sharepaxdist];
z=[x, herfCityPair, totpotential, totsinglepot];			
y=totenter;
w=inv(z'*z);

%Simulation draws are taken for each market and each airline in each market
simDraws=10;
t=simDraws+simDraws*numAirlines;
%sim=normrnd(0,1,rows(x),t);
%save sim sim
load sim;		%sim is a list of random variables as large as 30
sim=sim(:,1:t);


%To save computational time I arbitrarily picked a value for rho
rho=.5;


%Specify the data needed for the function
data.x=x;
data.xi=xi;
data.z=z;
data.sim=sim;
data.simdraws=simDraws;
data.numairlines=numAirlines;
data.rho=rho;
data.y=y;
data.w=w;

options = optimset('TolFun',.00000000001,'TolX',.00000001,'MaxIter',5000);


%BerryGMMsim2([-2;.2;.2;.7;10],data);


%I use the gridSearch function to get near the estimated value.  Once I'm near the estimated 
%value I can use the Nelder-Mead simplex search method as in Berry's paper.  The
%objective function is very nonlinear and discontinuous.
tic
f=gridSearch('BerryGMMsim2',[-.55,-.61;1.5,1.65;.4,.8;1.2498,1.2498;7.7448,7.7448],[5;5;15;1;1],data,0,0)
f.b
f.val
toc

%-0.5950    1.5375    0.5429    1.2498    7.7448
%  271.8204
[b,o,j,c]=repeatSearch('BerryGMMsim2',f.b,options,data,8,.00000000001);
b
o

 %-0.5963    1.5459    0.5547    1.2470    7.7243
 %270.3395
f.b=[-.5963 1.5459  0.5547 1.247  7.7243]'; 


%I calculat the variance covariance matrix twice to show that how the
%numerical derivative is taken may change the estimates of the standard
%errors.
var=numVarianceCov('BerryGMMsim2MOMENT',f.b,data,w,.01);
var.se

var=numVarianceCov('BerryGMMsim2MOMENT',f.b,data,w,.02);
var.se



