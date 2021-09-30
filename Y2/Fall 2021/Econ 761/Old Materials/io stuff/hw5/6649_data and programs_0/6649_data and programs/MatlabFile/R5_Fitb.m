function [yMean yCov Pi SmChurn]= R5_Fitb(par,const,E,Esm,data,nei_id, nei_dist)
% Report the fit of the model: 
%   model mean and sample mean, and corr bw model prediction and sample
% par=[pop_k,py_k,urb_k,mw_k,const_k,dkw,dk,dks;
%      pop_w,py_w,urb_w,dist,south,const_w,dwk,dw,dws,rho;
%      pop_s,py_s,urb_s,south_s,const_s,dsk,dsw,dss; tao, fc];
% const is a struct, with the following fields:
% c.N: sample size;
% c.NR: number of repetitions
% c.maxsm: max number of retailers

% data=[CTID,LnPop,LnRT,Urb,MW,   LnDist,South,KmS,WmS,Ns,  
%   Tk,Tw,KmS0WmS,NeiK,NeiW;    LnPop78,LnRT77,Urb80,Ns78]
% nei_id is the CTID for adjacent cnties; nei_dist is the inverse of
%     distance in miles

N=const.N; NR=const.NR; maxsm=const.maxsm;
oneC=ones(N,1); oneR=ones(1,NR); oneRS=ones(1,maxsm);
if size(par,1)==1
    par=par';
end;

%%% Create parameter vector and regressors
bk = par(1:5); dk=par(7);
bw = par(9:14); dw=par(16);
rho = par(18);            %0<rho<1
sig=realsqrt(1-rho^2);

bs=par(19:23); dss=par(26);
dpar=[par(6:8); par(15:17); par(24:26)]; 
tao=par(27);    %persistency of the market error
sc=par(28);     %sunk cost
bs2=par(29);

Tk=data(:,11); Tw=data(:,12);
Xk = [data(:,2:5),oneC]*bk + dk*Tk;    %add in the neighbor effect of counties not in sample
Xw = [data(:,[2:4,6:7]),oneC]*bw + dw*Tw;    %add in the neighbor effect of counties not in sample
Xs=[data(:,[2:4,7]),oneC]*bs;

Eps=E(:,1:NR);
Eps_prime=tao*Eps+sqrt(1-tao^2)*E(:,3*NR+1:4*NR);
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);

% Obtain the number of model predicted sm in 1978
Xs78=[data(:,[16:18,7]),oneC]*[bs(1:end-1);bs2]+sc;    %entrants in 78 pay the sunk cost
tmp = Xs78(:,oneR) + Eps; %avoid repmat to increase speed
Sm78=floor(exp(-tmp/dss));      % NxNR matrix; don't be confused w/ Ns78
Sm78=min(Sm78,maxsm);

% average profit for small stores in 1978
PiSm78=zeros(NR,1);
for i=1:NR
    tmp=Xs78+Eps(:,i);
    ts=Sm78(:,i);
    ts(ts==0)=1;
    tmp=tmp+dss*log(ts);
    tmp(Sm78(:,i)==0)=[];
    PiSm78(i)=mean(tmp);
end;
PiSm78=mean(PiSm78);

% Only new entrants incur fixed cost in 88
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));

In78=zeros(size(Esm));  %indicator of whether open in 78
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    [tmp Indx]=sort(tmp,2,'descend');
    In78(:,(i-1)*maxsm+1:i*maxsm)= (Indx<=tsm(:,oneRS));
    
    tmp = tmp + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
clear E Eps* Esm CtM SCM tmp 

%simulated average of Nk,Nw,Ns, Nk*Nei, Nw*Nei
const.dpar=dpar;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
Pi= R5_SPNGamePi(const,XKM,XWM,XSM,nei_id,nei_dist,Tk,Tw);
Pi=[Pi(1:2);PiSm78;Pi(3:end)];
SmChurn=R5_SPNGameChurn(const,XKM,XWM,XSM,nei_id,nei_dist,In78);

% report the sample mean and the model predicted mean
yMean=zeros(6,2);   %Km,Wm,Ns78,Ns,NeiK,NeiW
yMean(:,1)=mean(data(:,[8:9,19,10,14:15]))';
yMean(:,2)=[mean(Nkws(:,1:2))';mean(mean(Sm78,2));mean(Nkws(:,[4,8:9]))'];

% report the correlation between model and sample
yCov=zeros(6,1);
yCov(1)=corr(Nkws(:,1),data(:,8));
yCov(2)=corr(Nkws(:,2),data(:,9));
yCov(3)=corr(mean(Sm78,2),data(:,19));  %Ns78
yCov(4)=corr(Nkws(:,4),data(:,10));     %Ns
yCov(5)=corr(Nkws(:,8),data(:,14));     %NeiK
yCov(6)=corr(Nkws(:,9),data(:,15));     %NeiW


    