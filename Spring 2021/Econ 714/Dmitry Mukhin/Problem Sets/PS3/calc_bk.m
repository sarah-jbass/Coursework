function [rhoI,htauI] = calc_bk(rhoI0,tol,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
    diff = 1;
    iter = 1;
    maxiter = 1e6;

    %Calculate difference between rhoI and rhoI0
    while (diff>tol&&iter<maxiter)
        [rhoI,htauI] = calc_bkrho(rhoI0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
            psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
        diff = abs(rhoI - rhoI0);
        iter = iter+1;
        rhoI0 = rhoI;
    end
end

%Calculate rhoI
function [rhoI,htauI] = calc_bkrho(rhoI0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
    htauI = calc_bkhtauI(rhoI0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    rhoI = htauI(1:end-1)\htauI(2:end);
end

%Calculate htauI
function htauI = calc_bkhtauI(rhoI,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    [A,B] = LOM(rhoI,rhoa,rhog,rhoL,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    [Q, Lambda] = eig(A);
    iQ = inv(Q);
    if abs(Lambda(1))>1; sel = 1 ; else; sel = 2; end
    C = iQ*B;
    lm = diag(Lambda);
    lam = lm(sel);
    Theta = (-1/lam)*C(sel,:)*inv(eye(4) - (1/lam)*diag([rhoa rhog rhoL rhoI]));
    zs = [a g htauL];
    vs = [k c];
    htauI = Theta(4)^(-1) * (vs* iQ(sel,:)' - (zs*Theta(1,1:3)'));
end
