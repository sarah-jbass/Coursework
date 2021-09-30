function [impres,alph_post,var_post,alph_pr,v_pr] = calc_post_impres_BVAR(Y,X,ident,pr_sel,bands,nsim,ihor,mp,scl,IV,lambda)
if strcmp(pr_sel,'NIW')
    hyp_nconj = define_nconj_hyp(Y,X);
    [impres,alph_post,var_post,alph_pr,v_pr] = draw_impres_NIW(Y,X,hyp_nconj,lambda,ident,bands,nsim,ihor,mp,scl,IV);
else % assume flat prior
    [impres,alph_post,var_post] = draw_impres_flat(Y,X,ident,bands,nsim,ihor,mp,scl,IV);
    alph_pr = []; v_pr = [];
end
end

function [impres,alphpost,Vpost] = draw_impres_flat(y,x,ident,bands,nsim,ihor,mp,scl,IV)
% This estimates the impulse responses as identified from either cholesky 
% or by IV, under a flat prior. Draws betas from multivariate t distribution.
% 7/31: also returns posterior of alpha, mean and variance
[T,Nx] = size(x);
[~,Ny] = size(y);
b_flat = x\y;
Sig = (y-x*b_flat)'*(y-x*b_flat)/(T-Nx);
alphpost = b_flat(:);
Vpost = (1/(T-Nx - Ny - 1))*kron(Sig.*(T-Nx),inv(x'*x));
ir_flat =zeros(1+ihor,Ny,nsim);
% iD = inv(diag(Vpost));
% crr = iD*Vpost*iD;
draws = mvtrnd(Vpost,T-Nx - Ny - 1,nsim);
for i=1:nsim
    bdraw = reshape(diag(Vpost).^(1/2).*draws(i,:)'+alphpost,Nx,Ny);
    rdraw = y-x*bdraw;
    if strcmp(ident,'IV')
        biv = IV\rdraw;
        sz = scl/biv(mp);
        shock = biv*sz;
    else % Cholesky 
    lc = chol((rdraw'*rdraw)/(T-Nx))';
    shock = lc(:,mp);
    sz = 1/shock(mp);
    shock = shock*sz;
    end
    ir_flat(:,:,i) = draw_ir(bdraw(1:end-1,:),ihor,shock);
end
ir_flat = sort(ir_flat,3);
ind = round(nsim*bands);
impres = ir_flat(:,:,ind);
end

function impres = draw_ir(B,irhor,shock)
%Draw impulse response for SVAR
ny = size(B,2);
impres = NaN(irhor+1,length(shock));
impres(1,:) = shock;
tempx = zeros(1,size(B,1));
tempx(1:ny) = shock(:);
for tt = 1:irhor
    impres(tt+1,:) = tempx*B;
    tempx(ny+1:end) = tempx(1:end-ny);
    tempx(1:ny) = impres(tt+1,:);
end
end


function [impres,apost,V_nconj,alphprior,vpr] = draw_impres_NIW(y,x,hyp_nconj,lambda,ident,bands,nsim,ihor,mp,scl,IV)
[T,Nx] = size(x);
[~,Ny] = size(y);
p = (Nx-1)/Ny;
b_ols = x\y;
Sig = (y-x*b_ols)'*(y-x*b_ols)/(T-Nx); % ols stuff
% define priors
alphprior = 0*b_ols(:);
for i=1:Ny
    alphprior(Nx*(i-1)+i) = hyp_nconj.alphprior(i);
end
Aprior = reshape(alphprior,Nx,Ny);
Vprior = (lambda^2)*hyp_nconj.Vprior; %eye(Nx)*hyp_nconj.Vprior;
iVprior = inv(Vprior);
nuprior = hyp_nconj.nuprior;
Sprior = hyp_nconj.Sprior;
vpr = kron(Sprior,Vprior);
Vpost = inv(iVprior + x'*x);
Apost = Vpost*(iVprior*Aprior + x'*x*b_ols);
apost = Apost(:);
Spost = Sig*(T-Nx) + Sprior + b_ols'*(x'*x)*b_ols + Aprior'*iVprior*Aprior - Apost'*(iVprior + x'*x)*Apost;
nupost = nuprior+T;
iSpost = inv(Spost);
iSpost = (iSpost+iSpost')./2; % this just ensures that this variable is symmetric - sometimes off due to numerical error
ir_nconj =zeros(1+ihor,Ny,nsim);
V_nconj = (1/(nupost - Nx - 1))*kron(Spost,Vpost);
V_nconj = (V_nconj + V_nconj')/2; % ensures symmetry
for i=1:nsim % This sequentially draws sigmas and then betas|sigmas. This is equivalent to drawing beta from its marginal.
    sigdraw = iwishrnd(iSpost,nupost); % draw sigma
    bdraw = reshape(mvnrnd(apost,kron(sigdraw,Vpost)),Nx,Ny); % draw coefficients|sigma
    rdraw = y-x*bdraw;
    
    %lc = chol((rdraw'*rdraw)/(T-Nx))';
    %lc = chol(sigdraw)';
    
    if strcmp(ident,'IV')
        biv = IV\rdraw;
        sz = scl/biv(mp);
        shock = biv*sz;
    else % Cholesky
        lc = chol(sigdraw)';
        shock = lc(:,mp);
        sz = scl/shock(mp);
        shock = shock*sz;
    end
    
    ir_nconj(:,:,i) = draw_ir(bdraw(1:end-1,:),ihor,shock);
    
end
ir_nconj = sort(ir_nconj,3);
ind = round(nsim*bands);
impres = ir_nconj(:,:,ind);
end

function hyp_nconj = define_nconj_hyp(Y,X) % this prior is closer to GLP than M-A & R
[T,Ny] = size(Y);
[~,Nx] = size(X);
sigsq = var(Y);
arSig = 0*sigsq;
for k=1:Ny % calculate ar(p) Sigma for each variable
    Xtemp = X(:,[ (k:Ny:Nx-1)  Nx]);
    rtemp = (eye(T) - Xtemp*inv(Xtemp'*Xtemp)*Xtemp')*Y(:,k);
    arSig(k) = (rtemp'*rtemp)/(T-size(Xtemp,2)); % unbiased estimate
end
hyp_nconj.alphprior = ones(Ny,1); % Only need to define first lag, all other lags are zero
hyp_nconj.Vprior = zeros(Nx,1);%zeros(Ny*Nx,1); % This is the main diagonal, all off-diagonals will be zero
for k=1:Nx % Set prior - own lags, other lags, etc
if k<Nx
    l = floor(k/Ny)+1; % lag
    o = mod(k,Ny);
    if o==0
        o=Ny;
        l=l-1;
    end
    hyp_nconj.Vprior(k) = 1/(arSig(o)*l^2);
else
    hyp_nconj.Vprior(k) = 10/min(sigsq);%10*max(sigsq); % very flat prior on intercept
end
end
hyp_nconj.Vprior = diag(hyp_nconj.Vprior);
hyp_nconj.nuprior = Ny+2;
hyp_nconj.Sprior = diag(arSig);
end