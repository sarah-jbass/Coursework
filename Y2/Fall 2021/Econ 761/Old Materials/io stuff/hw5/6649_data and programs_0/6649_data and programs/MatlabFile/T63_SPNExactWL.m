function WL= T63_SPNExactWL(dw,XW,SimW,nei_id,nei_dist,aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the lower bound of Wm's best response        %%%%%%%%%%%%%%%%%%%
% SimW is the result from last iteration

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
WL=SimW;
