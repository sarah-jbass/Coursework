%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate a 3 stage model
% Set the max number of potential small entrants to 11
% Allow for a different constant in year 78 for small stores
% second stage
% June 2nd, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step one: create data matrix for year 88
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load R5_3stage88T XMat nei_id nei_dist data78 Weight2  %XMat nei_id nei_dist data78
% XMat=[CTID,LnPop,LnRT87,Urb90,MW,  LnDist,South,Kmm,Wmm,Ns,  Tk,Tw,...
%    (Kmm==0).*Wmm, NeiK, NeiW];
% data78={'statecd','cntycd','lnpop78','lnpy81','urb80','sm78','lnrt_pc77','smmax1
% ','smmax2','urbpop80','py81','rt77','pop78','pop80';}

% const.N: sample size
% const.NR: number of repititions
% const.MaxSm: max number of small retailers (number of potentail entrants)
% const.prm: prime numbers
% const.state: state of the random generator
const.N=size(XMat,1);
const.NR=150;
const.maxsm=11;
const.prm=primes(100);
const.state=5489;
[E Esm]=R5_Halton(const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step two: first stage estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Par0=[1.4033;2.2001;2.2721;0.5229;-24.5917;-0.3292;0.5874;-0.0138;
    1.3914;1.6759;2.4203;-1.4921;1.0472;-10.7131;-1.1126;1.24459;-0.0181;
    0.6782;1.5277;1.1503;-1.4074;0.9211;-9.7057;-0.9907;-0.9306;-2.3122;0.5706;-1.7815;-8.6205];
%%% This is the estimate using identity matrix as weights

opts = optimset('fminsearch');
opts.Display = 'iter';
opts.MaxFunEvals = 3000;
opts.MaxIter=2000;

%%% Parameter: bk, dkw, dkk, dks; bw, dwk,dww,dws; rho; bs, dsk,dsw,dss
lb = [-inf; -inf; -inf; -inf; -inf; -inf; 0.001; -.25;
    -inf; -inf; -inf; -inf; -inf; -inf; -inf; 0.001; -inf; 0.05;
    -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; 
    0.05; -inf; -inf];
ub = [+inf; +inf; +inf; +inf; +inf; -0.15; +inf; -0.001; ...
    +inf; +inf; +inf; +inf; +inf; +inf; -0.001; +inf; -0.001; 0.95;
    +inf; +inf; +inf; +inf; +inf; -0.001; -0.001; -0.001;
    0.95; -0.001; +inf];

tic
[Par2, fval2,exitflag2,output2] = fminsearchbnd(@R5_SPNFn2,Par0,lb,ub,opts,const,E,Esm,...
    [XMat,data78(:,[3,7,5,6])],nei_id, nei_dist,Weight2,'R5_88T14_para')
[yMean yCov Pi SmChurn]= R5_Fitb(Par2,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id, nei_dist);
SimSm= R5_CompSm(Par2,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id, nei_dist);

np=size(Par2,1);
out1=zeros(np+2,2);
out1(1:np,:)=[Par0,Par2];
out1(np+1:end,2)=[fval2*2065;exitflag2];

out2=zeros(27,2);
out2(1:6,:)=yMean;
out2(7:end-2,2)=[yCov;Pi;SimSm;SmChurn];

report=[out1;zeros(1,2);out2];
clear E Esm
save R5_3stage88T14 

