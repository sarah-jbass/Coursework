clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 710\Mikkel Solvsten\Problem Sets\PS3'

%Written by Sarah Bass, 2/12/21

%Clean data
reload = 0;
if reload
    x1 = readtable('AK91.csv');
    y = x1.lwage;
    ed = x1.educ;
    n = length(y);
    yob = zeros(n,9);
    sob = zeros(n,50);
    qob = zeros(n,3);
    for t = 1:n
        if x1.yob(t)>30 && x1.yob(t) <40
            yob(t,x1.yob(t)-30) = 1;
        end
        if x1.sob(t)>0 && x1.sob(t)<51
            sob(t,x1.sob(t)) = 1;
        end
        if x1.qob(t)>1 && x1.qob(t)<5
            qob(t,x1.qob(t)-1) = 1;
        end
    end
    save 'cleandata'
else
    load 'cleandata'
end

%Only use states for which we have data
ssob = sum(sob);
ssobi = ssob>0;
sob = sob(:,ssobi);

%Vars
x = [ed ones(n,1) yob sob];
z = [qob ones(n,1) yob sob];

%Beta and residual
beta = (x'*z*inv(z'*z)*z'*x)\(x'*z*inv(z'*z)*z'*y);
e = y - x*beta;

qzz = (z'*z)/n;
qxz = (x'*z)/n;
omega = 0*qzz;
for i=1:n
    omega = omega + z(i,:)'*z(i,:)*e(i)^2;
end
omega = omega/n;

minv = inv(qxz*inv(qzz)*qxz');
vhb = minv*(qxz*inv(qzz)*omega*inv(qzz)*qxz')*minv/n;
sqrt = sqrt(vhb(1,1));


