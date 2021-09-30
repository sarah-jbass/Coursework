function [E Esm]=R5_Halton(const)
% create Halton sequences
% const.N: sample size
% const.NR: number of repititions
% const.maxsm: max number of small retailers (number of potentail entrants)
% const.prm: prm numbers
% const.state: state of the random generator

N=const.N; 
NR=const.NR;
maxsm=const.maxsm;
prm=const.prm;

% Set the random seeds
rand('twister', const.state);

% Generate the random component for market error, Km, Wm
E=zeros(N,NR*4);
E(:,1:NR)=T30_halton(prm(1),NR,N,191);        %Eps
E(:,NR+1:2*NR)=T30_halton(prm(2),NR,N,397);   %EtaK
E(:,2*NR+1:3*NR)=T30_halton(prm(3),NR,N,587);   %EtaW
E(:,3*NR+1:4*NR)=T30_halton(prm(4),NR,N,937);   %Eps_prime

% Random component for the small retailers: shuffled halton
Hsm=zeros(N,maxsm);        %Halton sequence for small retailers
for i=1:maxsm
    Hsm(:,i)=T30_halton(prm(i+4),1,N,37);   
end;

Esm=zeros(N,NR*maxsm);
for i=1:NR
    for k=1:maxsm
        tid=randperm(N);
        Esm(:,(i-1)*maxsm+k)=Hsm(tid,k);
    end;
end;

