%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate a 3 stage model
% Set the max number of potential entrants to 11
% Allow for a different constant in year 78 for small stores
% second stage
% June 2nd, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step one: create data matrix for year 97
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load R5_3stage97T XMat nei_id nei_dist data78 Weight2  %XMat nei_id nei_dist data78
% XMat=[CTID,LnPop,LnRT87,Urb90,MW,  LnDist,South,Kmm,Wmm,Ns,  Tk,Tw,...
%    (Kmm==0).*Wmm, NeiK, NeiW];
% data78={'statecd','cntycd','lnpop78','lnpy81','urb80','sm78','lnrt_pc77','smmax1
% ','smmax2','urbpop80','py81','rt77','pop78','pop80';}

% const.N: sample size
% const.NR: number of repititions
% const.MaxSm: max number of small retailers (number of potentail entrants)
% const.prm: prm numbers
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
Par0=[1.4998;2.1557;1.3397;0.4009;-24.2813;-0.7312;0.6137;-0.027;
    2.0317;1.9896;1.6652;-1.0516;0.9245;-16.9369;-0.6814;0.85663;-0.0976;0.68322;
    1.6351;1.3801;-1.8835;1.1402;-11.6334;-0.445;-0.7862;-2.6727;
    0.31471;-2.3334;-9.6889];
%%% This is the estimate using identity matrix as weights

opts = optimset('fminsearch');
opts.Display = 'iter';
%opts.TolX = 1.e-6;
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
    [XMat,data78(:,[3,7,5,6])],nei_id, nei_dist,Weight2,'R5_97T17')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate standard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const.NR=300;
const.NPar=size(Par2,1);
const.NG=size(GMat,2);
[E Esm]=R5_Halton(const);

GMat= R5_SPNGb(Par2,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id, nei_dist);
CovFn2=R1_CovFn(GMat,nei_id,KMN);
disp('Eigenvalues of V at 2nd-stage theta, nonnegative then posi-semi-defi')
eig(CovFn2)

%%% Find out the Gredient of the moment functions
G=mean(GMat)';
DG=zeros(const.NG,const.NPar);
for i=1:const.NPar
    pp = Par2;
    inc = pp(i)*0.001;
    if abs(inc)<0.0005
        inc=0.0005*sign(pp(i));
    end;
    pp(i)=pp(i)+inc;
    [tmp,tmp1]=R5_SPNFn2(pp,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id,nei_dist,Weight2,'temp');
    DG(:,i)=(tmp1'-G)/inc;
end;

%%% Find out the standard deviation of para estimates
Avar = DG' * Weight2 * DG;
Avar = inv(Avar)/const.N;
T = Par2./sqrt(diag(Avar));

Avar2 = DG' * inv(CovFn2) * DG;
Avar2 = inv(Avar2)/const.N;
T2 = Par2./sqrt(diag(Avar2));

disp('2nd stage Para, T, and T2')
[Par2, T, T2]

clear E Esm
save R5_3stage97T17 
diary off
