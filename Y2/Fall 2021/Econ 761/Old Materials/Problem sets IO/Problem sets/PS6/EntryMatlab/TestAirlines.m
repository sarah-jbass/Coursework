

diary 'C:\Documents and Settings\Abe\Desktop\Entry\probit\results';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Examine entry in the airline industry using the probit model


%Load the data
clear all
run loadLongData
iota = ones(length(enter),1);
X=[iota pop city2 distance];

simDraws=200;
u=normrnd(0,1,rows(X),simDraws);
data.y=enter;
data.x=X;
data.z=X;
data.sim=u;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The first probit is adapted from Mike Cliff's GMM estimation
%Specify options for the GMM function below
op=foptions;
op(14)=8000;
op(2)=.000001;
op(3)=.000000000001;


b0 = zeros(size(X,2),1);
clear opt in
opt.gmmit = 1;
opt.W0 = 'I';
opt.plot = 0;
opt.vname = strvcat('Const','pop','city2','distance');

opt.infoz.momt = 'probitm';
opt.infoz.jake = 'probitj';
opt.S = 'probitS';
out3 = gmm(b0,opt,enter,X,X);


%Specify the options needed for the fminsearch program used below
op = optimset('TolFun',.00000000001,'TolX',.00000001,'MaxIter',5000);

%Specify my own GMM function and solve using the fminsearch program to 
%solve the GMM problem specified above.  Note that the estimation values 
%are the same.
[f1,o1]=fminsearch('probitGMMfunction',out3.b+[0.1;-0.2;0.6;0],op,data);
f1
01
%W=inv((1/size(X,1))*X'*X);
W=inv(X'*X);
%Note that the weight matrix here is not the efficient weight matrix.
%I have more work to do here, but the 
w=numVarianceCov('probitMomentFunc',f1,data,W,.0000001);
w.se


%Specify a simulation version of the probit model and use the fminsearch value
%to solve this problem near the solution.
[f3,o3,j]=repeatSearch('probitGMMsimulation',f1,op,data,20,.00000000001);
f3
o3
W=inv((1/size(X,1))*X'*X);
%Again, I need to do more work to solve for the standard errors.
w=numVarianceCov('probitSimulationMoment',f3,data,W,.01);
w.se



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search programs:
%gridSearch - Creates a grid and searches for points in the
%grid that minimize the objective function.
%
%repeatSearch - Repeatedly applies the Nelder-Mead search algorithm
%until the local minimum is reached.  The algorithm is repeatedly applied
%because it often gets stuck and must be restarted at the point it 
%last attempted.


tic
f4=gridSearch('probitGMMsimulation',[-3.5,-4;0,.5;3,4;.5,1],[5;5;5;5],data,0,0);
toc
f4.b											%Paramters on the grid that minimize the objective function
probitGMMsimulation(f4.b,data)		%Objective function at the value f4.b


[b4,o4,j,conv]=repeatSearch('probitGMMsimulation',f4.b,op,data,20,.00000000001);
b4												%Paramters after the repeatSearch
o4												%Value of objective function

%%%%%%
%Note that the above search still does not reach the value of the objective function
%when the true paramter is used.  
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section estimates the Berry model, but ignores the order of entry effects.
%Although we obtain a unique estimate of the number of firms, N, in the market, the 
%estimates do not pertain to a unique equilibrium.  It is shown in Berry's paper that
%an order of entry is imposed then the estimates pertain to a unique equilibrium solution.



%Load the data and specify market data x, firm data xi, instruments z, and 
%market entry information y

clear all
run loadWideData;
const=ones(size(pop,1),1);
x=[const pop distance];
xi=[sharepaxdist];
z=[x, xi];			%Remember that the number of insturements must be greater than num of paramters.
y=totenter;
w=inv(z'*z);

%Simulation draws are taken for each market and each airline in each market
simDraws=10;
t=simDraws+simDraws*numAirlines;
%sim=normrnd(0,1,rows(x),t);
load sim

rho=.5;


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

tic
f=gridSearch('BerryGMMsim2',[-.5,-.7;1.3,1.6;.08,.10;1.09,1.4;9,11],[5;5;5;5;5],data,0,0)
f.b
f.val
toc


[b,o,j,c]=repeatSearch('BerryGMMsim2',f.b,options,data,8,.00000000001);
b
o
% -0.6007    1.5605    0.1042    1.2989   10.0329
%  335.8096


%Again, I need to do more work to solve for the standard errors.
var=numVarianceCov('BerryGMMsim2',f.b,data,w,.01);
var.se

diary off;


