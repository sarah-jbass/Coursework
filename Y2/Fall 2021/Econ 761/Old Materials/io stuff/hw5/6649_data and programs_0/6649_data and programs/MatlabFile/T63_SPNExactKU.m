function KU= T63_SPNExactKU(dk,XK,SimK,nei_id,nei_dist,aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the upper bound of Km's best response        %%%%%%%%%%%%%%%%%%%
% SimK is the result from last iteration
% aa is the max number of iteration

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
