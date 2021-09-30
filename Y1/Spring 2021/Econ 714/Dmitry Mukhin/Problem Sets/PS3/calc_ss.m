function [Ybar,Cbar,Kbar,Lbar] = calc_ss(palpha,psigma,pphi,pGbar,pAbar,ptaubarL,ptaubarI,pbeta,pdelta)
syms Y C K L
f1 = Y == pAbar * K^palpha *L^(1-palpha);
f2 = Y == C + pdelta*K + pGbar*Y;
f3 = L^pphi * C^psigma == (1-ptaubarL)*pAbar*(1-palpha)*K^palpha*L^(- palpha);
f4 =(1+ptaubarI) == pbeta*(pAbar*palpha*K^(palpha-1)*L^(1-palpha) + (1-pdelta)*(1+ptaubarI));
soln = solve([f1 f2 f3 f4],[Y C K L]);

%Check if real and >0
validY = zeros(length(soln.Y),1);
validL = validY;
validC = validY;
for j=1:length(soln.Y)
    if (isreal(eval(soln.Y(j,1)))&& eval(soln.Y(j,1))>0) 
        validY(j)=1;
    end
    if (isreal(eval(soln.L(j,1)))&& eval(soln.L(j,1))>0)
        validL(j)=1;
    end
    if (isreal(eval(soln.C(j,1)))&& eval(soln.C(j,1))>0)
        validC(j)=1;
    end
end

ind = find((validY>0).*(validL>0).*(validC>0)); %Where everything is >0
Ybar = eval(soln.Y(ind(1),1));
Cbar = eval(soln.C(ind(1),1));
Kbar = eval(soln.K(ind(1),1));
Lbar = eval(soln.L(ind(1),1));
end

        