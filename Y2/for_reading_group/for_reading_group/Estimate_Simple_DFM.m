function [DF,iter,params,params0,diff] = Estimate_Simple_DFM(YY,tune,init)
% Estimates dynamic factor model
[T,N] = size(YY);
T=T-1; % -1 because ar(1) is imposed
% Initialize guess for F
if nargin<3
pcw = pca(zscore(YY));
F0 = zscore(YY)*pcw(:,1);
else 
F0 = init;% = zscore(YY)*pcw(:,1); % you would normally just initialize w/o the noise
end
if nargin<2
    tune = 1;
end
% pcw = pca(zscore(YY));
% F0 = zscore(YY)*pcw(:,1);
% F0 = zscore(F0);
% assume F follows an ar(1) structure, no intercept. Initialize params
% object
params0.PHI = F0(1:end-1,:)\F0(2:end,:);
Fhat0 = F0(1:end-1,:)*params0.PHI; % B.2
params0.Q = (F0(2:end,:) - Fhat0)'*(F0(2:end,:) - Fhat0)/T;% MLE estimate of Q
params0.LAMBDA = F0\YY;%YY\F0; % B.3
params0.SIGMA_nu1 = (YY - F0*params0.LAMBDA)'*(YY - F0*params0.LAMBDA)/T;
params0.SIGMA_nu2 = params0.Q;
params = params0;

iter = 1; % Begin while loop
MaxIter = 1e6;
dtol = 1e-10;
diff = 99;
S10 = 0;
P10 = 1e3;
while (iter<MaxIter)&&(diff>dtol)
    %% Kalman filter to calculate F|params
    %F = kalman_filter(YY,params); % Kalmon filter (E step) to find F
    [~,F] = KFS(YY',[],params.PHI',1,[],params.LAMBDA',params.Q,params.SIGMA_nu1,S10,P10);
    F=F';
    %F=zscore(F');
    %if F(1)<0; F=-1.*F; end
    %figure
    %plot(F)
    diff = sum(abs(F-F0),1); % check for convergence
    F = tune*F + F0*(1-tune);
    % Given F, update params
    params.PHI = F(1:end-1,:)\F(2:end,:);
    Fhat = F(1:end-1,:)*params.PHI; % B.2
    params.Q = (F(2:end,:) - Fhat)'*(F(2:end,:) - Fhat)/T;% MLE estimate of Q
    params.LAMBDA = F\YY;%YY\F; % B.3
    params.SIGMA_nu1 = (YY - F*params.LAMBDA)'*(YY - F*params.LAMBDA)/T;
    params.SIGMA_nu2 = params.Q;
    F0 = F;
    iter = iter+1;
end
DF = F;
end

function Fac = kalman_filter(YY,params)
[TT,NN] = size(YY);
% kalman filter
S = zeros(TT,1); % mean S_t|Y^t-1
S_upd = S; % mean S_t|Y^t
P = S;  % var S_t|Y^t-1
P_upd = P; % var S_t|Y^t
Yhat = zeros(TT,NN); % mean Y_t|Y^t-1
F = zeros(TT,NN,NN); % var Y_t|Y^t-1
S_smooth = S;

S(1) = 0; % alpha_1|0
P(1) = 1e3;
x1 = params.LAMBDA';
x2 = 1;
D0 = 0;
D1 = params.PHI;
A = D1;
D2 = 1;

for tt=1:TT % At = D1; Bt = D2; Dt = x1; Ct = 0; Zt = 0; Qt = params.SIGMA_nu1; Rt = params.SIGMA_nu2; SIGMA(tt) = PT|T-1; SIGMA(tt) = PT|T
% expected squared error of Y_tt
F(tt,:,:) = params.LAMBDA*P(tt)*params.LAMBDA' +params.SIGMA_nu1;
% update state equations
S_upd(tt) = S(tt) + P(tt)*x1*inv(squeeze(F(tt,:,:)))*(YY(tt,:) - x1*S(tt))';
P_upd(tt) = P(tt) - P(tt)*x1*inv(squeeze(F(tt,:,:)))*x1'*P(tt);
% predict variables in state equation for next period
if tt<TT
    S(tt+1) = D1*S_upd(tt);
    P(tt+1) = D1*P_upd(tt)*D1' + D2*params.SIGMA_nu2*D2';
    Yhat(tt+1,:) = x1*S(tt+1);
end
end
S_smooth(end) = S_upd(end);
% smooth
for tt=TT-1:-1:1
    S_smooth(tt) = S_upd(tt) + P_upd(tt)*D1*inv(P(tt+1))*(S_smooth(tt+1) - S(tt+1));
end
Fac=S_smooth;
end