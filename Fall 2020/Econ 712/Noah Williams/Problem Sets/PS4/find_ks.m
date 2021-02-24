function[ks,V,c,aprime,dist,pr,pw] = find_ks(kd0,assets,tol,maxiter,palpha,pbeta,pgamma,pdelta,Q,l,V0)
pr = (palpha)*(kd0)^(palpha - 1); % FOC of production function
pw = (1-palpha)*(kd0)^(palpha)*(1)^(-palpha); % FOC of production fn %kd0^palpha - pr*kd0; %zero profit condition

na = length(assets);

legal = true(na,na,2); %identify assets that are feasible
for nl = 1:2
   for nap = 1:na
       legal(nap,:,nl) = (pw*l(nl) + (1+pr-pdelta)*assets(nap) - assets')>0; %VERIFY
   end
end

V=0*V0; %preallocation, value doesn't matter but same size as V
aprime = V; %preallocation
aind = V; %preallocation
diff=999; %reset difference
iter=1; %reset iterations
while ((iter<maxiter)&&(diff>tol))
    %disp(['beginning iteration ' num2str(iter) ', diff = ' num2str(diff)])
    for nl = 1:2 %each labor endowment
    for ai = 1:na %and each asset holding
        Val = (pw*l(nl) + (1+pr-pdelta)*assets(ai) - assets(legal(ai,:,nl))' ).^(1-pgamma)./(1-pgamma) ...
            + pbeta*(V0(legal(ai,:,nl),1)*Q(nl,1) + V0(legal(ai,:,nl),2)*Q(nl,2));
        if isempty(Val); Val = -9999999999; end
        [V(ai,nl),ind] = max(Val); %fill in V
        aprime(ai,nl) = assets(ind); %fill in assets w/ index for max(Val)
        aind(ai,nl) = ind; %fill in index
    end
    end
    iter = iter+1;
    diff = sum(abs(V(:)-V0(:))); %slower rate of convergence than norm, functionally the same
    V0 = V;
end

% Stationary distribution
dist0 = ones(size(V));
dist0 = dist0./(sum(dist0(:))); %every weight is 1/2na
dist = dist0;
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    dist = 0*dist0;
    for nl = 1:2
       for ai = 1:na
           if dist0(ai,nl)>1e-14 % 
               targ = aind(ai,nl);%assets==aprime(ai,nl);
               dist(targ,1) = dist(targ,1)+dist0(ai,nl)*Q(nl,1);
               dist(targ,2) = dist(targ,2)+dist0(ai,nl)*Q(nl,2);
           end
       end
    end
    iter = iter+1;
    diff = sum(abs((dist0(:)-dist(:))));
    dist0 = dist;
end
marginal = sum(dist,2);
ks = sum(marginal.*assets');
c  = pw*repmat(l,na,1) + pr*repmat(assets',1,2) -aprime + (1-pdelta)*repmat(assets',1,2);
end

%Checks!
% sum(dist(:,1)) %should =.25
% sum(dist(:,2)) %should = .75
% sum(dist(:)) %should =1