clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 709\Jack Porter\Problem Sets\PS4'

%Written by Sarah Bass, 12/7/20

%Import Data
[data,sh] = xlsread('cps09mar.xlsx','Sheet1');
sh = sh(1,:);

%Define vars
n_y = 'earnings';
i_y = strcmp(n_y,sh);
i_h = strcmp('hours',sh);
i_w = strcmp('week',sh);
y = log(data(:,i_y)./(data(:,i_h).*data(:,i_w)));
t = length(y);
n_x = {'education' 'age'};
x=ones(t,1);
for i=1:length(n_x)
    x = [x data(:,strcmp(n_x{i},sh))];
end
x =[x (x(:,1+find(strcmp(n_x,'age'))) - x(:,1+find(strcmp(n_x,'education'))) -6)];
x(:,1+find(strcmp(n_x,'age'))) = x(:,end);
x(:,end) = x(:,end).^2/100;
n_x = {'education' 'experience' 'experience\^{}2/100' 'constant'};

%Single Asian Males
marital = data(:,strcmp(sh,'marital'));
race = data(:,strcmp(sh,'race'));
female = data(:,strcmp(sh,'female'));
male = ~female;
single = marital==7; %is this right?
asian = race==4;
sam = logical(single.*asian.*male);
x =x(sam,:);
y = y(sam,:);
weirdobs = x(:,3)>=45;
x = x(~weirdobs,:);
y = y(~weirdobs,:);

% reorder to put constant last
x = [x(:,2:end) x(:,1)];
beta = (x'*x)\(x'*y);

%Problem Set 4 material
r = y-x*beta;
[n,k] = size(x);
xpxi = inv(x'*x);

% compute robust standard errors
lev = 0*y;
for i=1:n
    lev(i) = x(i,:)*xpxi*x(i,:)';
end
xxpr = 0*xpxi;
lxxpr = xxpr;
lsxxpr = xxpr;
for i=1:n
    xxpr = xxpr + x(i,:)'*x(i,:)*r(i)^2;
    lxxpr = lxxpr + (1-lev(i))^(-1)*x(i,:)'*x(i,:)*r(i)^2;
    lsxxpr = lsxxpr + (1-lev(i))^(-2)*x(i,:)'*x(i,:)*r(i)^2;
end
HC0 = xpxi*xxpr*xpxi; 
se0 = sqrt(diag(HC0));
disp('SE0 estimates')
disp(num2str(se0'))

t=table([beta(1); se0(1)], [beta(2); se0(2)], [beta(3); se0(3)], [beta(4); se0(4)],'RowNames',{'Coefficient' 'Robust SE'},'VariableNames',{'Edu' 'Exp' 'Exp\^{}2/100' 'Constant'});
%table2latex(t,'table_hw4')
that = beta(1)/(beta(2) + beta(3)/5);
disp(['theta hat: ' num2str(that)] )

den = beta(1)/that;
gp = [1/den; -beta(1)/(den^2); -(beta(1)/5)/(den^2)];
sth = sqrt(gp'*HC0(1:end-1,1:end-1)*gp);
disp(['s(theta): ' num2str(sth)] )
disp(['CI: [' num2str(that - sth) ',' num2str(that + sth) ']'])