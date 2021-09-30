clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 714\Dmitry Mukhin\Problem Sets\PS3'

%Written by Sarah Bass, 2/9/21

%Prepare Data
clean=0;
if clean
    reload = 1;
    if reload
        [x,xt] = xlsread('macropset3q3.xls', 'matlab');
        [I,It] = xlsread('macropset3q3.xls', 'inv');
        save 'raw'
    else 
        load 'raw'
    end 
    names = xt(1,2:end); % First row of the xt matrix, stored as a cell (for strings)
    dates = datetime(xt(2:end,1)); %Formatting as dates
    dateI = datetime(It(2:end,1));
    
    %Question 1 - plot x
    figure
    for i= 1:3
        subplot(3,1,i)
        plot(dates,x(:,i),'k')
        title(names{i})
    end
    set(gcf,'Color',[1,1,1])
    saveas(gcf, 'rawdata.png')
    
    %Plot investment
    figure
    plot(dateI,I(:,1),'k')
    title('Investment')
    set(gcf,'Color',[1,1,1])
    saveas(gcf, 'rawi.png')
    
    %Question 2
    x = log(x);
    I = log(I);
    xhp = hpfilter(x); %Hodrick Prescott filter
    Ihp = hpfilter(I);

    figure
    for i=1:3
        subplot(3,1,i)
        hold on
        plot(dates,x(:,i),'k')
        plot(dates,xhp(:,i),'r')
        hold off
        title(names{i})
    end
    legend('log data','HP trend')
    set(gcf,'Color',[1 1 1])
    saveas(gcf,'log.png')

    figure
    hold on
    plot(dateI,I,'k');
    plot(dateI,Ihp,'r');
    hold off
    set(gcf,'Color',[1 1 1])
    title('Log investment data with trend')
    saveas(gcf,'logi.png')

    x = x-xhp;
    I = I - Ihp;

    figure
    for i=1:3
        subplot(3,1,i)
        plot(dates,x(:,i),'k')
        title(names{i})
        ylabel('log deviation from trend')
    end
    set(gcf,'Color',[1 1 1])
    saveas(gcf,'det.png')

    figure
    plot(dateI,I,'k');
    ylabel('log deviation from trend')
    set(gcf,'Color',[1 1 1])
    title('Log detrended investment data')
    saveas(gcf,'deti.png')

    %Question 3
    pdelta = 0.025;
    k = 0*I;
    for i = 2:length(k)
        k(i) = (1-pdelta)*k(i-1) + pdelta*I(i-1); %Linearized LOMK
    end
    T = size(x,1);
    k = k(end-T+1:end); %Only need capital from 1980 onwards

    %Time series of capital
    figure 
    plot(dates,k(:,1),'k')
    title('Capital Approximation')
    set(gcf,'Color',[1,1,1])
    saveas(gcf, 'lomk.png')

    close all 
    I = I(end-T+1:end);
    y = x(:,1);
    c = x(:,2);
    l = x(:,3);
    save 'clean'
    else 
        load 'clean'
end 

%Question 4
%Solve for steady states
palpha = 1/3;
psigma = 1;
pphi = 1;
pGbar = 1/3;
pAbar = 1;
ptaubarL = 0;
ptaubarI = 0;
pbeta = 0.99;

a = y - palpha*k - (1-palpha)*l;
[Ybar,Cbar,Kbar,Lbar] = calc_ss(palpha,psigma,pphi,pGbar,pAbar,ptaubarL,ptaubarI,pbeta,pdelta);
Ibar = pdelta*Kbar;
g = (1/pGbar)*(y - (Cbar/Ybar)*c - (Ibar/Ybar)*I);
htauL = (-1)*(pphi*l + psigma*c - palpha*k + palpha*l);

%Calculate shock persistance
rhoa = a(1:end-1)\a(2:end);
rhog = g(1:end-1)\g(2:end);
rhoL = htauL(1:end-1)\htauL(2:end);

%% Blanchard-Kahn Calculations
rhoI0 =0;
tol = 1e-6;
rerunbk = 1;
if rerunbk
    [rhoI,htauI] = calc_bk(rhoI0,tol,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,...
        psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    save 'runbk'
else
    load 'runbk'
end

figure
hold on
plot(dates,a,'k')
plot(dates,g,'r')
plot(dates,htauL,'b')
plot(dates,htauI,'m')
hold off
set(gcf,'Color',[1 1 1])
title('Wedges since 1980')
legend('a','g','\tau_L','\tau_I')
saveas(gcf,'wedges.png')

figure
subplot(2,2,1)
plot(dates,a,'k')
title('a')
subplot(2,2,2)
plot(dates,g,'k')
title('g')
subplot(2,2,3)
plot(dates,htauL,'k')
title('\tau_L')
subplot(2,2,4)
plot(dates,htauI,'k')
title('\tau_I')
set(gcf,'Color',[1 1 1])
saveas(gcf,'wedges2.png')

%% Counterfactual

reruncounter = 1;
if rerunbk
    [ga,gg,ghtauL,ghtauI] = counterfact_bk(rhoI,rhoa,rhog,rhoL,a,g,htauL,htauI,...
        c,k,palpha,pdelta,psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    save 'counterbk'
else
    load 'counterbk'
end

%Plots
figure
hold on
plot(dates,ga,'k--')
plot(dates,gg,'r--')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m--')
plot(dates,y,'Color',[0 0.6 0])
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
saveas(gcf,'wedgesg.png')

date2007 = dates(datenum(dates)>datenum('12/31/2007'));
inds = find(datenum(dates)>datenum('12/31/2007'));
date2007 = datetime(date2007(datenum(date2007)<datenum('12/31/2009')));
inds = inds(datenum(date2007)<datenum('12/31/2009'));

figure
hold on
plot(dates,ga,'k--')
plot(dates,gg,'r--')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m--')
plot(dates,y,'Color',[0 0.6 0])
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2007(1) date2007(end)])
saveas(gcf,'wedgesfin.png')

dga = ga(inds,:) - ga(inds(1),:);
dgg = gg(inds,:) - gg(inds(1),:);
dghtauL = ghtauL(inds,:) - ghtauL(inds(1),:);
dghtauI = ghtauI(inds,:) - ghtauI(inds(1),:);
dy = y(inds,:) - y(inds(1),:);

figure
hold on
plot(date2007,dga,'k--')
plot(date2007,dgg,'r--')
plot(date2007,dghtauL,'b--')
plot(date2007,dghtauI,'m--')
plot(date2007,dy,'Color',[0 0.6 0])
hold off
set(gcf,'Color',[1 1 1])
ylabel('Net change since 2008')
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2007(1) date2007(end)])
saveas(gcf,'wedgesfindiff.png')

date2020 = dates(datenum(dates)>datenum('12/31/2019'));
inds = find(datenum(dates)>datenum('12/31/2019'));
date2020 = datetime(date2020(datenum(date2020)<datenum('12/31/2020')));
inds = inds(datenum(date2020)<datenum('12/31/2020'));

figure
hold on
plot(dates,ga,'k--')
plot(dates,gg,'r--')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m--')
plot(dates,y,'Color',[0 0.6 0])
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2020(1) date2020(end)])
saveas(gcf,'wedgescov.png')

dga = ga(inds,:) - ga(inds(1),:);
dgg = gg(inds,:) - gg(inds(1),:);
dghtauL = ghtauL(inds,:) - ghtauL(inds(1),:);
dghtauI = ghtauI(inds,:) - ghtauI(inds(1),:);
dy = y(inds,:) - y(inds(1),:);

figure
hold on
plot(date2020,dga,'k--')
plot(date2020,dgg,'r--')
plot(date2020,dghtauL,'b--')
plot(date2020,dghtauI,'m--')
plot(date2020,dy,'Color',[0 0.6 0])
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
ylabel('Net change since 2020')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2020(1) date2020(end)])
saveas(gcf,'wedgescovdiff.png')
