function [A,B] = LOM(rhoI,rhoa,rhog,rhoL,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    syms C K Cp Kp I Ip L Lp g a htauL htauI gp ap htauLp htauIp
    
    %Functions
    f1 = ap == rhoa*a;
    f2 = gp == rhog*g;
    f3 = htauLp == rhoL*htauL;
    f4 = htauIp == rhoI*htauI;
    f5 = Kp == (1-pdelta)*K + pdelta*I;
    f6 = a + palpha*K + (1-palpha)*L == (Cbar/Ybar)*C + (pdelta*Kbar/Ybar)*I + pGbar*g;   
    f7 = ap + palpha*Kp + (1-palpha)*Lp == (Cbar/Ybar)*Cp + (pdelta*Kbar/Ybar)*Ip + pGbar*gp;
    f8 = pphi*L + psigma*C == (-1)*htauL + palpha*K - palpha *L;
    f9 = pphi*Lp + psigma*Cp == (-1)*htauL + palpha*Kp - palpha *Lp;
    f10 = psigma*(Cp - C) + htauI == pbeta*(palpha*pAbar*Kbar^(palpha - 1)*Lbar^(1-palpha)*(ap + (1-palpha)*(Lp - Kp)) +(1-pdelta)*htauIp);
    
    %Solve
    obj = solve([f1 f2 f3 f4 f5 f6 f7 f8 f9 f10], [Cp Kp I Ip L Lp gp ap htauLp htauIp]);
    exp = [obj.Kp; obj.Cp];
    X = [K C];
    Z = [a g htauL htauI];
    A = eval(jacobian(exp,X));
    B = eval(jacobian(exp,Z));
end