function f=probitMomentFunc(b,data);


y=data.y;
x=data.x;
z=data.z;
%function [m,e] = probitm(b,infoz,stat,y,x,z)
%
% Moments forProbit Model.  See gmmldv_d for a demo.
  
q = 2*y-1;			% Convert from 0/1 to -1/+1 coding
F = normcdf(x*b.*q);
f = normpdf(x*b.*q);
f.e = q.*f./F;



f.N=length(z);
f.moment=(1/f.N)*z'*f.e;

L=size(z,2);
temp=(f.e*ones(1,L)).*z;
f.omega=(1/f.N).*temp'*temp;


