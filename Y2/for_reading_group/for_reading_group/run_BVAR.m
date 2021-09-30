clear; close all; clc
% This file generates the BVAR figures from my presentation.
% Most of this file is for plotting, most of the actual BVAR estimation is
% in calc_post_impres_BVAR.m
% Note: In these codes I am assuming more observations than variables in X
% To do otherwise you need to impose the prior using dummy observations
% This is straightforward but not how I wrote these codes
load 'vardata'
m = {'m0' 'm1' 'm2' 'm3' 'm4' 'm5' 'm6' 'm7' 'm8' 'm9' 'm10'};
for im = 1:length(m)
    if strcmp(m{im},'m1') % base
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';% 'IV', 'Cholesky'
        excl_covid = 1; 
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 1;
    elseif strcmp(m{im},'m0') % base, but flat prior (just for likelihood graph not impulse response graph)
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';% 'IV', 'Cholesky'
        excl_covid = 1; 
        pr_sel = 'flat'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 1;
    elseif strcmp(m{im},'m2') % base, closed economy
        var_incl = {'US IP' 'US RS' 'US CPI' ...
            'US M2' 'US Credit (% GDP)' 'FFR'  ...
            'US HY Spread' '10 yr Treasury' ...
            'Effective Bond Premium' }; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 3 0 0 0 0];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 8;
        lambda = 1;
    elseif strcmp(m{im},'m3') % base
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 0.1;
    elseif strcmp(m{im},'m4') % base, closed economy
        var_incl = {'US IP' 'US RS' 'US CPI' ...
            'US M2' 'US Credit (% GDP)' 'FFR'  ...
            'US HY Spread' '10 yr Treasury' ...
            'Effective Bond Premium' }; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 3 0 0 0 0];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 8;
        lambda = 0.1;
    elseif strcmp(m{im},'m5') % base
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 1000;
    elseif strcmp(m{im},'m6') % base, closed economy
        var_incl = {'US IP' 'US RS' 'US CPI' ...
            'US M2' 'US Credit (% GDP)' 'FFR'  ...
            'US HY Spread' '10 yr Treasury' ...
            'Effective Bond Premium' }; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 3 0 0 0 0];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 8;
        lambda = 1000;
    elseif strcmp(m{im},'m7') % flat prior
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'flat'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 1;
    elseif strcmp(m{im},'m8') % flat prior, closed
        var_incl = {'US IP' 'US RS' 'US CPI' ...
            'US M2' 'US Credit (% GDP)' 'FFR'  ...
            'US HY Spread' '10 yr Treasury' ...
            'Effective Bond Premium' }; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 3 0 0 0 0];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'IV';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'flat'; % 'NIW' 'flat'
        p = 6;
        mp = 8;
        lambda = 1;
    elseif strcmp(m{im},'m9') % cholesky
        var_incl = {'US IP' 'US RS' 'UK IP' 'DE IP' 'US CPI' ...
            'US M2' 'UK CPI' 'US Credit (% GDP)' 'FFR' 'UK Interbank' 'DE Interbank'  ...
            'US HY Spread' '10 yr Treasury' 'Real RiskAversion' 'GlobRiskAsPr' ...
            'Effective Bond Premium' 'Dollars/Euro' 'Dollars/Pound'}; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 2 2 2 0 0 0 0 0 0 0 3 0 2 2];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'Cholesky';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 13;
        lambda = 1;
    elseif strcmp(m{im},'m10') % cholesky, closed
        var_incl = {'US IP' 'US RS' 'US CPI' ...
            'US M2' 'US Credit (% GDP)' 'FFR'  ...
            'US HY Spread' '10 yr Treasury' ...
            'Effective Bond Premium' }; % order matters for Cholesky identification, but not IV
        tr = [ 2 2 2 2 3 0 0 0 0];% 0 is demean, 1 is log demean, 2 is log detrend
        ident = 'Cholesky';%'Cholesky'; % 'IV', 'Cholesky'
        excl_covid = 1; % cut sample off at end-2019
        pr_sel = 'NIW'; % 'NIW' 'flat'
        p = 6;
        mp = 8;
        lambda = 1;
    end
    
    [T,~] = size(raw);
    nvars = length(var_incl);
    clean = zeros(T,nvars);
    Xtemp = [ones(T,1) (1:T)']; 
    for k=1:nvars % Transform and include the selected variables
        temp = raw(:,strcmp(var_incl{k},Varnames));
        switch tr(k)
            case 0
                temp = temp - mean(temp);
            case 1
                temp = 100*(log(temp) - mean(log(temp)));
            case 2
                temp = 100*((eye(T) - Xtemp*inv(Xtemp'*Xtemp)*Xtemp')*log(temp));
            case 3
                temp = (eye(T) - Xtemp*inv(Xtemp'*Xtemp)*Xtemp')*temp;
        end
        clean(:,k) = temp;
    end
    if excl_covid % exclude covid?
        clean = clean(year(dates)<2020,:);
        cleandates = dates(year(dates)<2020,:);
        IV = MPshock(year(dates)<2020,:);
    else
        cleandates = dates;
        IV = MPshock;
    end
    if im == 2 || im ==3 % plots
    figure('Renderer', 'painters', 'Position', [10 10 1110 860])
    for k=1:nvars
        subplot(ceil(nvars/3),3,k)
        plot(cleandates,clean(:,k),'k')
        title(var_incl{k})
    end
    set(gcf,'Color',[1 1 1])
    suptitle('Variables entering VAR')
    cd('pings')
    saveas(gcf,['model_' num2str(im) '.png'])
    cd('..')
    end
    % Make X and Y matrices
    X = lagmatrix(clean,1:p);
    X = X(p+1:end,:);
    Y = clean(p+1:end,:);
    IV = IV(p+1:end,:);
    [T,Nx] = size(X);
    X = [X ones(T,1)]; T=T+1;
    bands = [0.1 0.5 0.9]; % 80% Coverage Intervals
    nsim = 1000; % number of draws for bands
    ihor = 24; % IR horizon
    scl = 0.25; % 25 basis point shock to 'mp' variable
    [impres,alph_post,var_post,alph_pr,v_pr] = calc_post_impres_BVAR(Y,X,ident,pr_sel,bands,nsim,ihor,mp,scl,IV,lambda); 
    
    xx = 0:ihor;
    if ~strcmp(m{im},'m0') % more plots
    figure('Renderer', 'painters', 'Position', [10 10 1110 860])
    for k=1:nvars
    subplot(ceil(nvars/3),3,k)
    plot(xx,impres(:,k,2),'k')
    hold on
    plot(xx,impres(:,k,1),'b-.')
    plot(xx,impres(:,k,3),'b-.')
    plot(xx,0*xx,'k:')
    hold off
    if k==1
        ylabel('%')
    end
    title(var_incl{k})
    end
    set(gcf,'Color',[1 1 1])
    legend('Median','80% CI','Location','SouthEast')
    if strcmp(pr_sel,'NIW')
    suptitle(['Shock identified by ' ident ', prior: ' pr_sel ', \lambda = ' num2str(lambda)])
    else
    suptitle(['Shock identified by ' ident ', prior: ' pr_sel])
    end
    cd('pings')
    saveas(gcf,['VAR_' m{im} '.png'])
    cd('..')
    else
        flat_alph = alph_post;
        flat_V = var_post;
    end
    % plot prior/posterior distributions
    if (strcmp(m{im},'m1')||strcmp(m{im},'m3'))||strcmp(m{im},'m5')
        xx = -0.5:0.01:1.5;
        figure('Renderer', 'painters', 'Position', [10 10 1110 860])
        subplot(4,1,1) % first lag
        plot(xx,normpdf(xx,alph_pr(1),sqrt(v_pr(1))),'k') % prior
        hold on
        plot(xx,normpdf(xx,alph_post(1),sqrt(var_post(1))),'r-.') % post
        plot(xx,normpdf(xx,flat_alph(1),sqrt(flat_V(1))),'b--') % flat post
        hold off
        ylim([0 18])
        legend('Prior','Posterior','Flat posterior')
        title('First lag')
        subplot(4,1,2) % second lag
        plot(xx,normpdf(xx,alph_pr(1+nvars),sqrt(v_pr(1+nvars,1+nvars))),'k') % prior
        hold on
        plot(xx,normpdf(xx,alph_post(1+nvars),sqrt(var_post(1+nvars,1+nvars))),'r-.') % post
        plot(xx,normpdf(xx,flat_alph(1+nvars),sqrt(flat_V(1+nvars,1+nvars))),'b--') % flat post
        hold off
        ylim([0 18])
        title('Second lag')
        subplot(4,1,3) % third lag
        plot(xx,normpdf(xx,alph_pr(1+2*nvars),sqrt(v_pr(1+2*nvars,1+2*nvars))),'k') % prior
        hold on
        plot(xx,normpdf(xx,alph_post(1+2*nvars),sqrt(var_post(1+2*nvars,1+2*nvars))),'r-.') % post
        plot(xx,normpdf(xx,flat_alph(1+2*nvars),sqrt(flat_V(1+2*nvars,1+2*nvars))),'b--') % post
        hold off
        ylim([0 18])
        title('Third lag')
        subplot(4,1,4) % third lag
        plot(xx,normpdf(xx,alph_pr(1+3*nvars),sqrt(v_pr(1+3*nvars,1+3*nvars))),'k') % prior
        hold on
        plot(xx,normpdf(xx,alph_post(1+3*nvars),sqrt(var_post(1+3*nvars,1+3*nvars))),'r-.') % post
        plot(xx,normpdf(xx,flat_alph(1+3*nvars),sqrt(flat_V(1+3*nvars,1+3*nvars))),'b--') % post
        hold off
        ylim([0 18])
        title('Fourth lag')
        set(gcf,'Color',[1 1 1])
        suptitle(['US Industrial Production: Prior and posterior: \lambda = ' num2str(lambda)])
        cd('pings')
        saveas(gcf,['pr_post_' m{im} '.png'])
        cd('..')
    end
end