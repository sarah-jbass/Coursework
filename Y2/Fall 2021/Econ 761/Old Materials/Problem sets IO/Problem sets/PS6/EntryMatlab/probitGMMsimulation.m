
%%%%
%This function calculates the moment function of a simulated
%probit model.
%%%%

function m= probitGMMsimulation(b,data)

y=data.y;
x=data.x;
z=data.z;
sim=data.sim;

%%%
%Calculate the predicted values of the function given the simulation draws.
%%%

temp=(x*b(1:length(b))*ones(1,size(sim,2))+sim)>=0;


%Average the simulation draws and compute how this differs
%from the observed value.

e=y-sum(temp,2)/size(sim,2);


%Calculate the value of the moment function.
m=(e'*z)*eye(size(z,2))*(z'*e);
