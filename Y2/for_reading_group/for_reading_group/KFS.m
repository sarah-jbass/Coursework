function [LL,S_sm] = KFS(Y,Z,A,B,C,D,Q,R,S10,P10)
% Kalman filter/smoother
% Nattinger, 7/2021
% Given Y,Z data and initializations for S1|0, P1|0, and params ABCDQR,
% Calculates posterior estimates of S, given by the smoothed S_sm
% Also calculates the log likelihood of the system, LL
[Ny,T] = size(Y); % figure out sizes for preallocation, and preallocate
Ns = size(A,1);
S_sm = zeros(Ns,T);
S_pr = S_sm;
S_upd = S_pr;
P_pr = zeros(Ns,Ns,T);
P_upd = P_pr;
F = zeros(Ny,Ny,T);
Y_pr = zeros(Ny,T);

S_pr(:,1) = S10; % initializations
P_pr(:,:,1) = P10;
for t=1:T % iterate through the system, one time index at a time
    if isempty(Z)
        Y_pr(:,t) = D*S_pr(:,t);
    else
        Y_pr(:,t) = C*Z(:,t) + D*S_pr(:,t);
    end
    F(:,:,t) = D*reshape(P_pr(:,:,t),Ns,Ns)*D' + R;
    S_upd(:,t) = S_pr(:,t) + reshape(P_pr(:,:,t),Ns,Ns)*D'*inv(reshape(F(:,:,t),Ny,Ny))*(Y(:,t) - Y_pr(:,t));
    P_upd(:,:,t) = reshape(P_pr(:,:,t),Ns,Ns)*D'*inv(reshape(F(:,:,t),Ny,Ny))*D*reshape(P_pr(:,:,t),Ns,Ns);
    
    if t<T
        S_pr(:,t+1) = A*S_upd(:,t);
        P_pr(:,:,t+1) = A*reshape(P_upd(:,:,t),Ns,Ns)*A' + B*Q*B';
    end
end
S_sm(:,end) = S_upd(:,end); % initialize smoothed state
for t=T-1:(-1):1 % work backwards to smooth, one t at a time
    S_sm(:,t) = S_upd(:,t) + reshape(P_upd(:,:,t),Ns,Ns)*A'*inv(reshape(P_pr(:,:,t+1),Ns,Ns))*(S_sm(:,t+1) - S_upd(:,t+1));
end
LL=0; % calculate log likelihood of the system
for t=1:T
LL =  LL+ (log(det(reshape(F(:,:,t),Ny,Ny)) +(Y(:,t) - Y_pr(:,t))'*inv(reshape(F(:,:,t),Ny,Ny))*(Y(:,t) - Y_pr(:,t))));
end
LL = (-1/2)*LL - T*Ny*log(2*pi)/2;
end