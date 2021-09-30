function [KBR,Bgp] = T63_SPNKBR(dk,XK,KU,nei_id,nei_dist,aa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Km's exact best response        %%%%%%%%%%%%%%%%%%%
% KU is the result from the previous iteration, an upper bound

Bgp=[];

N=size(KU,1);           %Need to update KU because XK has changed
SimK=KU;
k=1;
while k<aa;
    tmpK=SimK;
    tmp=[SimK;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tk=tmp(nei_id);       
    Pi=XK+2*dk*sum(Tk.*nei_dist,2);
    
    SimK= (Pi>=0);    
    if tmpK==SimK
        break
    end;
    k=k+1;
end;
KU=SimK;

SimK=zeros(N,1);
k=1;
while k<aa;
    tmpK=SimK;
    tmp=[SimK;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tk=tmp(nei_id);       
    Pi=XK+2*dk*sum(Tk.*nei_dist,2);
    
    SimK= (Pi>=0);    
    if tmpK==SimK
        break
    end;
    k=k+1;
end;
KBR=SimK;               %The lower bound of true Km's Best Response, needs to be improved

%%% Step one: For those with KBR=0 and KU=1, check for pairs
T=find(KBR<KU);     %T are the obs with KBR=0 and KU=1
NT = size(T,1); 

if NT<2     %NT should only be 0. It guards against the case of NT=1 (which should not happen)
    return      %return to the program that calls T63_SPNKBR
else
    t_Pi=Pi(T);
    for i=1:aa
        K0=KBR(T);
        %%% Check for pairs first
        [KBR(T),t_Pi]=T63_SPNPair(T,nei_id(T,:),nei_dist(T,:),t_Pi,dk,aa);
        if K0==KBR(T)
            break
        end;

        T=find(KBR<KU);
        if numel(T)<2
            break
        end;
    end;
%     if i==aa
%         disp('pair iteration did not converge for KBR')
%     end;
end;

%%% Step two: for those with KBR=0 and KU=1, split into 'connected' groups
T=find(KBR<KU);
NT1=size(T,1);      

if NT1<=2       %if NT1=2, then the two obs with KBR=0 and KU=1 should remain 0
    return;     
else
    %%% Update Pi after Step One
    tmp=[KBR;0];      %the biggest index in nid is 2066, which indicates a empty value
    Tk=tmp(nei_id);       
    Pi=XK+2*dk*sum(Tk.*nei_dist,2);
    Ngp=T63_SPNNei(T,nei_id(T,:));
    
    t_N=size(Ngp,1);      %the number of elements in the connected group
    for i=1:t_N
        T=Ngp{i};
        if numel(T)>=3 && numel(T)<=7      %I have omitted the very few groups bigger than 7; 
                                           %in other applications, should increase the value 
           %any connected group with less than 2 obs should remain 0
            [MaxVec,maxtpi] =T63_SPNPi(T,nei_id(T,:),nei_dist(T,:),Pi(T),dk);
            if maxtpi>0
                KBR(MaxVec)=1;           %KBR is the optimal solution
            end;
        elseif numel(T)>7
%            disp('KBR conn grp >7')
            Bgp=[Bgp;T];
        end;
    end;
end;




