function  g=BerryGMMsim2(db,data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform the data into the different components
%y = the entry data
%x = the market level data
%xi = firm data
%z = instruments used for the gmm simulation
y=data.y;
x=data.x;
xi=data.xi;
z=data.z;
sim=data.sim;
simdraws=data.simdraws;
numairlines=data.numairlines;
rho=data.rho;
w=data.w;
N=length(y);

%Seperate the parameter values
%delta = parameter capturing the competitive effect
%bmkt = parameters for market level data
%bind = parameters for individual firm level data
delta=db(1,:);
bmkt=db(2:(1+size(x,2)));
bind=db(size(x,2)+2:size(db));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Profits:

%market profits for each firm = market variables + competitive effects + simulation random errors
mkt_profits=kron(x*bmkt,ones(1,simdraws))+rho*sim(:,1:simdraws);

%Expand the vector of individual firm coefficients so they are as large as the
%number of firms.  Multiply individual firm coefficients times corresponding
%firm data.  The data in tempMult is arranged as
%[b1.*(first factor of all airlines),b2.*(second factor of all airlines)] 
bindall=kron(bind,ones(numairlines,1));
tempMult=xi.*(ones(N,1)*bindall');

%For each airline sum up the components of individual profits.
%temp contains the estimated individual portion 
estindprofit=zeros(N,numairlines);
for i=1:length(bind)
	   estindprofit=estindprofit+tempMult(:,numairlines*(i-1)+1:numairlines*i);
end

%Calculate the individual portion of profit by adding the error term.
indpart_profits=kron(ones(1,simdraws),estindprofit)+(1-rho^2)^.5*sim(:,simdraws+1:simdraws*(numairlines+1));
%Calculate the total profits by adding the individual profits plus the market profits
totind_profits=kron(mkt_profits,ones(1,numairlines))+indpart_profits;


%Estimate the entry decision in each market given the simulation draws.
%For each simulation draw and for each market, we calculate whether
%it is feasible for a certain potential number of firms to enter profitably.
%Need to fix this so the potential number of entrants is not always 42.
pot=max(y);
entry=zeros(N,simdraws*pot);
k=vectorCombinations([pot;simdraws]);
temp=zeros(N,1);
for i=1:simdraws*pot
   	temp=sum(totind_profits(:,numairlines*(k(i,2)-1)+1:numairlines*k(i,2))>=-delta*(k(i,1)),2)>=k(i,1);
   	entry(:,i)=temp*k(i,1);  
end

%For each simulation draw and for each market, we calculate the maximum
%number of firms that could profitably enter.
estentry=zeros(N,simdraws);
for j=1:simdraws
   estentry(:,j)=max(entry(:,(j-1)*pot+1:pot*j),[],2);
end


%Calculate the market level error by subtracing the actual entry level 
%from the predicted level of entry.
e=y-sum(estentry,2)./simdraws;


%Calculate the moment conditions and moment function.
moment=e'*z;
g=diag(moment*w*moment');

