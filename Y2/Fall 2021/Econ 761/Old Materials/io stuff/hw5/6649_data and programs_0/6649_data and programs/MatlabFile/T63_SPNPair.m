function [Sim,Pi]=T63_SPNPair(T,nei_id,nei_dist,XB,d,nr)
%%% check whether some pairs can be changed to 1 simultaneously

ND=size(nei_id,2); big=max([max(nei_id(:)),max(T)]);
pi=-ones(big,1)*999; pi(T)=XB;
pitmp=pi(nei_id);
pitmp=pitmp+XB(:,ones(1,ND))+2*d*nei_dist;   %change in pi if entering market (i,j)
[t1, tmp]=find(pitmp>=0);     %t1 are the obs that should be changed to 1
t1=T(t1);
t2=nei_id(pitmp>=0);          %t2 are the neighbors that should be changed to 1
t=union(t1,t2); t(t==0)=[]; t=t(:); %all obs in t should be changed to 1

%%% Update the above pair-result
Sim=zeros(big,1);
Sim(t)=1;
k=1;
while k<nr
    tmp=Sim(T);
    Tn=Sim(nei_id);
    Pi=XB+2*d*sum(Tn.*nei_dist,2);
    
    Sim(T)=(Pi>=0);
    if tmp==Sim(T)
        break
    end;
    k=k+1;
end;

Sim=Sim(T);
Pi(Pi>=0)=[];