diary R5_3stage88T14Res.out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate variance of para in R5_3stage88T14
% June 18, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short g

load R5_3stage97T17 XMat
XMat97=XMat;

load R5_3stage88T14 Par2 XMat data78 nei_id nei_dist Weight2 const
const.N=size(XMat,1);
const.NR=300;
const.maxsm=11;
const.prm=primes(100);
const.state=5489;
[E Esm]=R5_Halton(const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate standard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const.NPar=size(Par2,1);
KMN=.5;

GMat= R5_SPNGb(Par2,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id, nei_dist);
CovFn2=R1_CovFn(GMat,nei_id,KMN);
disp('Eigenvalues of V at 2nd-stage theta, nonnegative then posi-semi-defi')
[min(eig(CovFn2)) max(eig(CovFn2))] 

%%% Find out the Gredient of the moment functions
G=mean(GMat)';
const.NG=size(GMat,2);
DG=zeros(const.NG,const.NPar);
for i=1:const.NPar
    pp = Par2;
    inc = pp(i)*0.001;
    if abs(inc)<0.0002
        inc=0.0002*sign(pp(i));
    end;
    
    pp(i)=pp(i)+inc;
    [tmp,tmp1]=R5_SPNFn2(pp,const,E,Esm,[XMat,data78(:,[3,7,5,6])],nei_id,nei_dist,Weight2,'temp');
    DG(:,i)=(tmp1'-G)/inc;
end;

%%% Find out the standard deviation of para estimates
Avar = DG' * Weight2 * DG;
Avar = inv(Avar)/const.N;
StdP=sqrt(diag(Avar));
T = Par2./StdP;

Avar2 = DG' * inv(CovFn2) * DG;
Avar2 = inv(Avar2)/const.N;
StdP2=sqrt(diag(Avar2));
T2 = Par2./StdP2;

disp('2nd stage Para, T, and T2')
StdP=[Par2, T, T2]

save R5_3stage88T14Res 
diary off
