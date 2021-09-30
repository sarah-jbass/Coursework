function [DF,DF_regions,iter,params,params0,diff] = Estimate_Regional_DFM(YY,tune,regions)
% Estimates dynamic factor model
[T,N] = size(YY);
T=T-1; % -1 because ar(1) is imposed
% Initialize guess for F
pcw = pca(zscore(YY));
F0 = zscore(YY)*pcw(:,1);
% Partition YY
ur = unique(regions);
nr = length(ur);
F_regions0 = zeros(T+1,nr);
YY0 = YY;
YY = [];
nir = zeros(1,nr);
for i=1:nr
    nir(i) = length(find(regions==ur(i)));
    Y_region = YY0(:,regions==ur(i));
    YY = [YY Y_region];
    Y_r_R = (eye(T+1) - F0*inv(F0'*F0)*F0')*Y_region;
    pcw = pca(zscore(Y_r_R));
    F_regions0(:,i) = zscore(Y_r_R)*pcw(:,1); % initial guess of regional factor is pca of Y in region, residualized w.r.t global factor
end
F0 = [F0 F_regions0];
if sum(nir,2)<N; error('Regions improperly defined'); end


% assume F follows a VAR(1) structure, no intercept. Initialize params
% object
params0.PHI = (F0(1:end-1,:)\F0(2:end,:))';
Fhat0 = F0(1:end-1,:)*params0.PHI; % B.2
params0.Q = (F0(2:end,:) - Fhat0)'*(F0(2:end,:) - Fhat0)/T;% MLE estimate of Q
%% This part changes due to zero restrictions
params0.LAMBDA = zeros(N,nr+1);
cnt = 1;
for i=1:nr
%params0.LAMBDA = YY\F0; % B.3
params0.LAMBDA(cnt:cnt+nir(i)-1,[1 1+i]) = ((F0(:,[1 1+i]))\YY(:,cnt:cnt+nir(i)-1))';%YY(:,cnt:cnt+nir(i)-1)\(F0(:,[1 1+i])); % restricted regressions
cnt = cnt+nir(i);
end

params0.SIGMA_nu1 = (YY - F0*params0.LAMBDA')'*(YY - F0*params0.LAMBDA')/T;
params0.SIGMA_nu2 = params0.Q;
params = params0;

iter = 1; % Begin while loop
MaxIter = 1e6;
dtol = 1e-1;
diff = 99;
S10 =[0;0;0];
P10 = eye(3)*1000;
while (iter<MaxIter)&&(diff>dtol)
    %% Kalman filter to calculate F|params
    %F = kalman_filter(YY,params,nir); % Kalmon filter (E step) to find F
    [~,F] = KFS(YY',[],params.PHI,1,[],params.LAMBDA,params.Q,params.SIGMA_nu1,S10,P10);
    F=F';
    diff = sum(sum(abs(F-F0),1),2); % check for convergence
    F = tune*F + F0*(1-tune);
    % Given F, update params
    params.PHI = (F(1:end-1,:)\F(2:end,:))';
    %
    params.LAMBDA = zeros(N,nr+1);
    cnt = 1;
    for i=1:nr
    params.LAMBDA(cnt:cnt+nir(i)-1,[1 1+i]) = YY(:,cnt:cnt+nir(i)-1)\(F(:,[1 1+i])); % restricted regressions
    cnt = cnt+nir(i);
    end
    Fhat = F0(1:end-1,:)*params.PHI; % B.2
    params.Q = (F(2:end,:) - Fhat)'*(F(2:end,:) - Fhat)/T;% MLE estimate of Q
    params.SIGMA_nu1 = (YY - F*params.LAMBDA')'*(YY - F*params.LAMBDA')/T;
    params.SIGMA_nu2 = params.Q;
    F0 = F;
    iter = iter+1;
end
DF = F(:,1);
DF_regions = F(:,2:end);
end

function Fac = kalman_filter(YY,params,nir)
[TT,NN] = size(YY);
nr = length(nir);
% kalman filter
S = zeros(TT,1+nr); % mean S_t|Y^t-1
S_upd = S; % mean S_t|Y^t
P = zeros(TT,1+nr,1+nr);  % var S_t|Y^t-1
P_upd = P; % var S_t|Y^t
Yhat = zeros(TT,NN); % mean Y_t|Y^t-1
F = zeros(TT,NN,NN); % var Y_t|Y^t-1
S_smooth = S;

S(1,:) = [0 0 0]; % alpha_1|0
P(1,:,:) = eye(3)*1e3;
x1 = params.LAMBDA';
x2 = 1;
D0 = 0;
D1 = params.PHI;
A = D1;
D2 = eye(3);

for tt=1:TT % At = D1; Bt = D2; Dt = x1; Ct = 0; Zt = 0; Qt = params.SIGMA_nu1; Rt = params.SIGMA_nu2; SIGMA(tt) = PT|T-1; SIGMA(tt) = PT|T
% expected squared error of Y_tt
F(tt,:,:) = params.LAMBDA*squeeze(P(tt,:,:))*params.LAMBDA' +params.SIGMA_nu1;
% update state equations
S_upd(tt,:) = S(tt,:) + (squeeze(P(tt,:,:))*x1*inv(squeeze(F(tt,:,:)))*(YY(tt,:) - S(tt,:)*x1)')';
P_upd(tt,:,:) = squeeze(P(tt,:,:)) - squeeze(P(tt,:,:))*x1*inv(squeeze(F(tt,:,:)))*x1'*squeeze(P(tt,:,:));
% predict variables in state equation for next period
if tt<TT
    S(tt+1,:) = (D1*S_upd(tt,:)')';
    P(tt+1,:,:) = D1*squeeze(P_upd(tt,:,:))*D1' + D2*params.SIGMA_nu2*D2';
    Yhat(tt+1,:) = S(tt+1,:)*x1;
end
end
S_smooth(end,:) = S_upd(end,:);
% smooth
for tt=TT-1:-1:1
    S_smooth(tt,:) = S_upd(tt,:) + (squeeze(P_upd(tt,:,:))*D1'*inv(squeeze(P(tt+1,:,:)))*(S_smooth(tt+1,:) - S(tt+1,:))')';
end
Fac=S_smooth;
end