function [WBR, Bgp]= T63_SPNWBR(dw,XW,WBR,nei_id,nei_dist,aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Wm's best response        %%%%%%%%%%%%%%%%%%%
% WBR is the result from the previous iteration, a lower bound

Bgp=[];

N=size(WBR,1);
SimW=ones(N,1);
k=1;
while k<aa;
    tmpW=SimW;
    tmp=[SimW;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tw=tmp(nei_id);       
    Pi=XW+2*dw*sum(Tw.*nei_dist,2);
    
    SimW= (Pi>=0);    
    if tmpW==SimW
        break
    end;
    k=k+1;
end;
WU=SimW;

SimW=WBR;
k=1;
while k<aa;
    tmpW=SimW;
    tmp=[SimW;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tw=tmp(nei_id);       
    Pi=XW+2*dw*sum(Tw.*nei_dist,2);
    
    SimW= (Pi>=0);    
    if tmpW==SimW
        break
    end;
    k=k+1;
end;
WBR=SimW;

%%% Step one: For those with WBR=0 and WU=1, check for pairs
T=find(WBR<WU);     %T are the obs with WBR=0 and WU=1
NT = size(T,1); 

if NT<2     %NT should only be 0. It guards against the case of NT=1 (which should not happen)
    return      %return to the program that calls T63_SPNWBR
else
    t_Pi=Pi(T);
    for i=1:aa
        W0=WBR(T);
        %%% Check for pairs first
        [WBR(T),t_Pi]=T63_SPNPair(T,nei_id(T,:),nei_dist(T,:),t_Pi,dw,aa);
        if W0==WBR(T)
            break
        end;

        T=find(WBR<WU);
        if numel(T)<2
            break
        end;
    end;
%     if i==aa
%         disp('pair iteration did not converge for WBR')
%     end;
end;

%%% Step two: for those with WBR=0 and WU=1, split into conn groups
T=find(WBR<WU);
NT1=size(T,1);      

if NT1<=2
    return;     
else
    %%% Update Pi after Step One
    tmp=[WBR;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tw=tmp(nei_id);       
    Pi=XW+2*dw*sum(Tw.*nei_dist,2);
    Ngp=T63_SPNNei(T,nei_id(T,:));
    
    t_N=size(Ngp,1);      %the number of connected group
    for i=1:t_N
        T=Ngp{i};
        if numel(T)>=3 && numel(T)<=7  %I omitted the very few groups bigger than 7
           %any connected group with less than 2 obs should remain 0
           [MaxVec,maxtpi] =T63_SPNPi(T,nei_id(T,:),nei_dist(T,:),Pi(T),dw);
            if maxtpi>0
                WBR(MaxVec)=1;           %WBR is the optimal solution
            end;
        elseif numel(T)>7
%            disp('WBR conn grp>7')
            Bgp=[Bgp;T];
        end;
    end;
end;

