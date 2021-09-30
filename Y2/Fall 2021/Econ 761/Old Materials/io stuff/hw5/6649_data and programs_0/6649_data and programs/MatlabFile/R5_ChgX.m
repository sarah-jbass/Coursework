function [ChgK ChgW ChgS ChgS78]=R5_ChgX(par,const,E,Esm,data,nei_id, nei_dist)
% Estimate the effect of changes in X
% par=[pop_k,py_k,urb_k,mw_k,const_k,dkw,dk,dks;
%      pop_w,py_w,urb_w,dist,south,const_w,dwk,dw,dws,rho;
%      pop_s,py_s,urb_s,south_s,const_s,dsk,dsw,dss; tao, fc,const_s78];
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
const.dpar=[par(6:8); par(15:17); par(24:26)]; 
tao=par(27);    %persistency of the market error
sc=par(28);     %sunk cost
bs2=par(29);

Eps=E(:,1:NR);
Eps_prime=tao*Eps+sqrt(1-tao^2)*E(:,3*NR+1:4*NR);

ChgK=zeros(6,1);
ChgW=zeros(7,1);
ChgS78=zeros(7,1);
ChgS=zeros(8,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Kmart's X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xw = [data(:,[2:4,6:7]),oneC]*bw + dw*data(:,12);    %add in the neighbor effect of counties not in sample
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);

% Obtain the number of model predicted sm in 1978
Xs78=data(:,[16:18,7])*bs(1:end-1)+bs2+sc;    %entrants in 78 pay the sunk cost; a different constant
tmp = Xs78(:,oneR) + Eps; %avoid repmat to increase speed
Sm78=floor(exp(-tmp/dss));      % NxNR matrix; don't be confused w/ Ns78
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,[2:4,7]),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;

% Base case
Xk = [data(:,2:5),oneC]*bk + dk*data(:,11);    
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(1)=sum(Nkws(:,1));    %total no of Km
ChgW(1)=sum(Nkws(:,2));     %total no of Wm
ChgS(1)=sum(Nkws(:,4));     %total no of sm
ChgS78(1)=sum(mean(Sm78,2));

% Pop increase by 10%
Xk = [log(1.1)+data(:,2),data(:,3:5),oneC]*bk + dk*data(:,11);   
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(2)=sum(Nkws(:,1));    %total no of Km

% Retail Sales increase by 10%
Xk = [data(:,2),log(1.1)+data(:,3),data(:,4:5),oneC]*bk + dk*data(:,11);  
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(3)=sum(Nkws(:,1));    %total no of Km

% Urban ratio increase by 10%
Xk = [data(:,2:3),1.1*data(:,4),data(:,5),oneC]*bk + dk*data(:,11);   
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(4)=sum(Nkws(:,1));    %total no of Km

% MidWest=0
Xk = [data(:,2:4),zeros(N,1),oneC]*bk + dk*data(:,11);   
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(5)=sum(Nkws(:,1));    %total no of Km

% MidWest=1
Xk = [data(:,2:4),oneC,oneC]*bk + dk*data(:,11);   
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgK(6)=sum(Nkws(:,1));    %total no of Km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Walmart's X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xk = [data(:,2:5),oneC]*bk + dk*data(:,11);   
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);

% Pop increase by 10%
Xw = [log(1.1)+data(:,2),data(:,[3:4,6:7]),oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(2)=sum(Nkws(:,2));    %total no of Wm

% Retail Sales increase by 10%
Xw = [data(:,2),log(1.1)+data(:,3),data(:,[4,6:7]),oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(3)=sum(Nkws(:,2));    %total no of Wm

% Urban increase by 10%
Xw = [data(:,2:3),1.1*data(:,4),data(:,6:7),oneC]*bw + dw*data(:,12); 
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(4)=sum(Nkws(:,2));    %total no of Wm

% Distance increase by 10%
Xw = [data(:,2:4),log(1.1)+data(:,6),data(:,7),oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(5)=sum(Nkws(:,2));    %total no of Wm

% South=0
Xw = [data(:,[2:4,6]),zeros(N,1),oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(6)=sum(Nkws(:,2));    %total no of Wm

% South=1
Xw = [data(:,[2:4,6]),oneC,oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgW(7)=sum(Nkws(:,2));    %total no of Wm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Sm's X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xw = [data(:,[2:4,6:7]),oneC]*bw + dw*data(:,12);   
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);

% Pop increases by 10%
Xs78=[log(1.1)+data(:,16),data(:,[17:18,7])]*bs(1:end-1)+bs2+sc; 
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss)); 
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[log(1.1)+data(:,2),data(:,[3:4,7]),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(2)=sum(mean(Sm78,2));
ChgS(2)=sum(Nkws(:,4));    

% Retail Sales increases by 10%
Xs78=[data(:,16),log(1.1)+data(:,17),data(:,[18,7])]*bs(1:end-1)+bs2+sc;    
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss));  
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,2),log(1.1)+data(:,3),data(:,[4,7]),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(3)=sum(mean(Sm78,2));
ChgS(3)=sum(Nkws(:,4));    

% Urban increases by 10%
Xs78=[data(:,16:17),1.1*data(:,18),data(:,7)]*bs(1:end-1)+bs2+sc;   
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss)); 
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,2:3),1.1*data(:,4),data(:,7),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(4)=sum(mean(Sm78,2));
ChgS(4)=sum(Nkws(:,4));    

% South=0
Xs78=[data(:,16:18),zeros(N,1)]*bs(1:end-1)+bs2+sc;  
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss));
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,2:4),zeros(N,1),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(5)=sum(mean(Sm78,2));
ChgS(5)=sum(Nkws(:,4));    

% South=1
Xs78=[data(:,16:18),oneC]*bs(1:end-1)+bs2+sc;   
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss));
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,2:4),oneC,oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(6)=sum(mean(Sm78,2));
ChgS(6)=sum(Nkws(:,4));    

% fixed cost increases by 10%
Xs78=data(:,[16:18,7])*bs(1:end-1)+bs2+sc*1.1;   
tmp = Xs78(:,oneR) + Eps; 
Sm78=floor(exp(-tmp/dss)); 
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
Xs=[data(:,[2:4,7]),oneC]*bs;
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
XSM=zeros(size(Esm));
for i=1:NR
    tsm=Sm78(:,i);
    SCM=1.1*sc*( CtM > tsm(:,oneRS) );  %the first Sm78 stores do not incur sunk cost
    tmp=Xs + sig*Eps_prime(:,i);
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm) + SCM;
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS78(7)=sum(mean(Sm78,2));
ChgS(7)=sum(Nkws(:,4));    

% How many small stores if all pays sc in the second stage
Xs=[data(:,[2:4,7]),oneC]*bs;
XSM=zeros(size(Esm));
CtM=1:maxsm;
CtM=CtM(oneC,:);       %CtM: a count matrix, with 1-11 for each row
for i=1:NR
    tmp=Xs + sig*Eps_prime(:,i) + sc;
    tmp = tmp(:,oneRS) + rho*Esm(:,(i-1)*maxsm+1:i*maxsm);
    tmp = sort(tmp,2,'descend') + dss*log(CtM);   %if tmp+dsk*Km+dsw*Wm>0 then enter
    XSM(:,(i-1)*maxsm+1:i*maxsm) =tmp;
end;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     
ChgS(8)=sum(Nkws(:,4));    



