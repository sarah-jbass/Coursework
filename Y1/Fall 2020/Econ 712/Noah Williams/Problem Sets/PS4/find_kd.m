function [kd] = find_kd(kd0,assets,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V0,tol2,maxiter2,tune)
iter = 1;
diff = 999;
kd = kd0;
while ((iter<maxiter2)&&(diff>tol2))
    ks = find_ks(kd,assets,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V0);
    iter = iter+1;
    diff = abs(kd-ks);
    kd = tune*kd + (1-tune)*ks;
    if kd<1e-10; kd = 1; diff = 999; end
end