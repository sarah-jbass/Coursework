function Ngp=T63_SPNNei(T,nei_id)
%%% Find out which observations are connected neighbors
%%% T is the observations. For ex: T=[2,4,5,6,7]
%%% nei_id are T's neighbors, a matrix whose row number is the size of T
%%% IVID(T) returns the index number of T: IVID(T)=[1,2,3,4,5]

nei_id(ismember(nei_id,T)==0)=0;

big=max(T);
IVID=zeros(big,1); 
IVID(T)=1; IVID=cumsum(IVID); 

t=T(1);
T(1)=[];
N=size(T,1);

k=1;
for i=1:N
    tid=IVID(t);    %the index of obs t
    tn=nei_id(tid,:);
    tn=unique(tn);
    tn=tn(:);
    tn_gp=intersect(tn,T);
    if numel(tn_gp)>0
        t=[t;tn_gp];
        T=setdiff(T, tn_gp);
        if numel(T)==0
            t=sort(t);
            Ngp{k}=t;
            break;
        end;
    elseif numel(tn_gp)==0
        t=sort(t);
        Ngp{k}=t;
        k=k+1;
        t=T(1);
        T(1)=[];
        if numel(T)==0
            Ngp{k}=t;
            break
        end;
    end;
end;
    


