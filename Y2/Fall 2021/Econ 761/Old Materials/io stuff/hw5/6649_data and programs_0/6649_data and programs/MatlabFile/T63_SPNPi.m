function [MaxVec,max_tpi] =T63_SPNPi(T,nei_id,nei_dist,X,d)
%%% Find out the exact vector that delivers the max profit
%%% T is the group of counties, w/ 3-7 obs 
%%% X is the XB part of the profit vector
%%% d is the size of the network effect
%%% nei_id and nei_dist are their neigh's ID and dist;

N=size(T,1);
ti=1:1:N;

NeiM=[];
for i=3:N
    t=nchoosek(ti,i);
    t(:,i+1:N)=0;
    NeiM=[NeiM;t];      %a matrix with all the combinations
end;

NP=size(NeiM,1);
tpi=zeros(NP,1);        %profit for each combination
for i=1:NP
    ti=NeiM(i,:);       %the observations' indices
    ti(ti==0)=[];       %remove '0' from the index list
    t=T(ti);            %the true observations
    tn=nei_id(ti,:);
    td=nei_dist(ti,:);
    td(ismember(tn,t)==0)=0;    %only change the observatios in t to 1
    tp=X(ti)+d*sum(td,2);
    tpi(i)=sum(tp);
end;
[max_tpi,max_id]=max(tpi); %MaxPi is the maximum profit; MaxTP is the corresponding vector
max_t=NeiM(max_id,:);
max_t(max_t==0)=[];
MaxVec=T(max_t);           %the combination of obs that achieves highest profit 

