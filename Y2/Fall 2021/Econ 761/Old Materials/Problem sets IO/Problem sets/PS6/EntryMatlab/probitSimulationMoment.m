function f=probitSimulationMoment(b,data);


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

f.e=y-sum(temp,2)/size(sim,2);

f.N=size(z,1);
f.moment=(1/f.N)*z'*f.e;



L=size(z,2);
temp=f.e*(ones(1,L)).*z;
f.omega=(1/f.N).*temp'*temp;


