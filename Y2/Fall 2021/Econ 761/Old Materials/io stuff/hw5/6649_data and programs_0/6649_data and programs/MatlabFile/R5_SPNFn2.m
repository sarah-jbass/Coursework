function [Moment G]= R5_SPNFn2(par,const,E,Esm,data,nei_id, nei_dist, Weight,outF)
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

Xk = [data(:,2:5),oneC]*bk + dk*data(:,11);    %add in the neighbor effect of counties not in sample
Xw = [data(:,[2:4,6:7]),oneC]*bw + dw*data(:,12);    %add in the neighbor effect of counties not in sample
Xs=[data(:,[2:4,7]),oneC]*bs;

Eps=E(:,1:NR);
Eps_prime=tao*Eps+sqrt(1-tao^2)*E(:,3*NR+1:4*NR);
XKM = Xk(:,oneR) + sig*Eps_prime + rho*E(:,NR+1:2*NR);
XWM = Xw(:,oneR) + sig*Eps_prime + rho*E(:,2*NR+1:3*NR);

% Obtain the number of model predicted sm in 1978
Xs78=data(:,[16:18,7])*bs(1:end-1)+bs2+sc;    %entrants in 78 pay the sunk cost; 
tmp = Xs78(:,oneR) + Eps; %avoid repmat to increase speed
Sm78=floor(exp(-tmp/dss));      % NxNR matrix; don't be confused w/ Ns78
Sm78=min(Sm78,maxsm);

% Only new entrants incur fixed cost in 88
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
clear E Eps* Esm CtM SCM tmp 

%simulated average of Nk,Nw,Ns, Nk*Nei, Nw*Nei
const.dpar=dpar;
Nkws= R5_SPNGame(const,XKM,XWM,XSM,nei_id, nei_dist,Sm78);     

%%% Define the moment functions
KmS=data(:,8); WmS=data(:,9); Ns=data(:,10); 
KmS0WmS=data(:,13); NeiK=data(:,14); NeiW=data(:,15);
Ns78=data(:,19);

Resid=[KmS,WmS,KmS0WmS, Ns, Ns.*KmS, Ns.*WmS, Ns.*KmS0WmS, NeiK, NeiW, Ns.*Ns78]-Nkws; 
Resid(:,11)= Ns78 - mean(Sm78,2);
Resid(:,12)= Resid(:,4)-Resid(:,11);    %dif in Ns between two periods
difX=data(:,2:3)-data(:,16:17);

data23=data(:,2).*data(:,3); data24=data(:,2).*data(:,4);
G=[Resid(:,1)'*[oneC,data(:,2:5),data23],...
    Resid(:,2)'*[oneC,data(:,[2:4,6:7]),data23],...
    Resid(:,3)'*[oneC,data(:,2:3),data23,data24],...
    Resid(:,4)'*[oneC,data(:,[2:4,7]),data23],...
    Resid(:,5)'*[oneC,data(:,2:3)],...
    Resid(:,6)'*[oneC,data(:,2:3)],...
    Resid(:,7)'*[oneC,data23,data24],...
    Resid(:,8)'*[oneC,data(:,2)],...
    Resid(:,9)'*[oneC,data(:,2)],...
    Resid(:,10)'*[oneC,data(:,[2:4,7])],...   %Ns78 .* Ns88 interacts w/ lnpop,lnrt,urb
    Resid(:,11)'*[oneC,data(:,16:18)],...   %Ns78 interacts w/ lnpop78,lnrt77,urb80
    Resid(:,12)'*difX ];                    %Ns-Ns78 interacts w/ lnpop-lnpop78,lnrt-lnrt77
G=G/N;

Moment = G*Weight*G';


%%% Save the iteration results
persistent iterCnt iterOut
 
if isempty(iterCnt)
    iterCnt=1;
else
    iterCnt=iterCnt+1;
end

if mod(iterCnt, 30)==1
    iterOut=[iterOut;
        iterCnt,par',Moment];
    eval(['save ' outF '.mat iterOut'])
end;

