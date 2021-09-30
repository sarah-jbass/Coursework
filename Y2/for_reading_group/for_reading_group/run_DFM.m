clear; close all; clc
% The following codes produce the DFM figures used in my presentation
% It looks like a lot but most of the code in this file are for plotting
% The DFM estimation occurs in: Estimate_Simple_DFM.m;
% Estimate_Regional_DFM.m.
load 'facdata'
T = size(XX,1);
cd('pings')
figure
plot(plotdates,zscore(XX))
set(gcf,'Color',[1 1 1])
title('Risky Asset Prices (standardized)')
ylabel('Z score')
saveas(gcf,'raw_ap.png')
YY = zscore(XX(2:end,:) - XX(1:end-1,:));
pcw = pca(YY);

figure
plot(plotdates(2:end),zscore((XX(2:end,:) - XX(1:end-1,:))))
set(gcf,'Color',[1 1 1])
title('Risky Asset Price Changes (standardized)')
ylabel('Z score')
saveas(gcf,'raw_ap_ch.png')


PCA_nosum = zscore((YY*pcw(:,1)));
figure
bar(plotdates(2:end),zscore(PCA_nosum),'r')
set(gcf,'Color',[1 1 1])
title('Risky Asset Price Changes PCA (standardized)')
ylabel('Z score')
saveas(gcf,'pca_ap_changes.png')

PCA = zscore(cumsum(YY*pcw(:,1)));
figure
plot(plotdates,zscore(mean(zscore(XX),2)),'k')
hold on
plot(plotdates(2:end),PCA,'r')
hold off
legend('Mean of standardized data (restandardized)','PCA factor (restandardized)','Location','SouthEast')
set(gcf,'Color',[1 1 1])
title('Aggregated Risky Asset Prices (standardized)')
ylabel('Z score')
saveas(gcf,'pca_ap.png')
cd('..')

% Estimate simple DFM
init = randn(T-1,1)*2 + YY*pcw(:,1); % to make things interesting, initialize with noise so that initial factor isn't correct

[DF,iter,params,params0,diff] = Estimate_Simple_DFM(YY,1,init);


cd('pings')
figure
plot(plotdates(2:end),zscore(PCA),'r')
hold on
plot(plotdates(2:end),zscore(cumsum(init)),'k')
plot(plotdates(2:end),zscore(cumsum(DF)),'b')
hold off
legend('PCA','Initialization (PCA+noise)','Dynamic Factor','Location','SouthEast')
set(gcf,'Color',[1 1 1])
title('Aggregated Risky Asset Prices (standardized)')
ylabel('Z score')
saveas(gcf,'DFM_ap.png')
cd('..')

% More interesting DFM, with regional factors in EM/AEs
EMs = ones(1,size(XX,2));
EMs([5 6 7 8 9 11 15 17 18 19 21 25 26 27 28 30 31 32 34 36 37 39 40 41 43]) = 0;

[DF_r,DF_regions,iter,params,params0,diff] = Estimate_Regional_DFM(YY,0.1,EMs);

cd('pings')
figure
plot(plotdates(2:end),zscore(PCA),'r')
hold on
plot(plotdates(2:end),zscore(cumsum(DF)),'b')
plot(plotdates(2:end),zscore(cumsum(DF_r)),'k')
hold off
legend('PCA','DFM, no regions','DFM, regions','Location','SouthEast')
set(gcf,'Color',[1 1 1])
title('Aggregated Risky Asset Prices (standardized)')
ylabel('Z score')
saveas(gcf,'DFM_comp.png')
cd('..')


%% Residualize off realized variance
GlobRiskAsPr = zscore([0; cumsum(DF_r)]); % Global risky asset prices - global factor estimated with regional factors
XX = [ones(T,1) RealVol];
RiskAversion = -1*zscore((eye(T) - XX*inv(XX'*XX)*XX')*GlobRiskAsPr);
figure
subplot(3,1,1)
plot(plotdates,GlobRiskAsPr,'k')
title('Global Risky Asset Prices')
subplot(3,1,2)
plot(plotdates,RiskAversion,'k')
title('Global Risk Aversion')
subplot(3,1,3)
plot(plotdates,RealVol*100^2,'k')
title('Realized Volatility')
set(gcf,'Color',[1 1 1])
suptitle('Asset Prices, Risk Aversion, and Volatility')
cd('pings')
saveas(gcf,'risk_aversion.png')
cd('..')

%% Misspecification? Real risky asset prices and risk aversion
xx=log(CPI(2:end,:)./CPI(1:end-1,:));
real_rap_ch = (eye(T-1) - xx*inv(xx'*xx)*xx')*DF_r; % deflate in changes using OLS - no intercept because we want to keep signs of inflation/deflation always consistent
real_rap = [0; cumsum(real_rap_ch)]; % convert to levels
real_RA_2 = -1*zscore((eye(T) - XX*inv(XX'*XX)*XX')*real_rap);
XX = [ones(T,1) RealVol CPI];
real_RiskAversion = -1*zscore((eye(T) - XX*inv(XX'*XX)*XX')*GlobRiskAsPr); % deflate in levels using OLS

figure
subplot(3,1,1)
plot(plotdates,RiskAversion,'k')
title('Nominal Risk Aversion')
subplot(3,1,2)
plot(plotdates,real_RiskAversion,'k')
title('Real Risk Aversion, AsPr deflated in levels')
subplot(3,1,3)
plot(plotdates,real_RA_2,'k')
title('Real Risk Aversion, AsPr deflated in changes')
set(gcf,'Color',[1 1 1])
suptitle('Nominal Risk Aversion vs Deflated Risk Aversion')
cd('pings')
saveas(gcf,'real_risk_aversion.png')
cd('..')