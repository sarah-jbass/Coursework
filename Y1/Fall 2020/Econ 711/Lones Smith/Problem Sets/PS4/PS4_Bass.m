clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 711\Lones Smith\Problem Sets\PS4'

%Written by Sarah Bass, 11/18/20

d = 0.05:0.05:10;
count=0;

for alpha = [0.25 0.5 0.75] %testing different values of alpha
        count=count+1;
        
    %Parameters
    aa = d./(2-alpha);
    ba = aa;
    apa = aa.*(d - aa + alpha*ba );%d.^2/(2-a)^2;
    bpa = ba.*(d - ba + alpha*aa );%APRA;
    tpa = apa + bpa;

    ab = d*(2+alpha)./(4 - 2*alpha^2);
    bb = (d./2)*(1+alpha*(2+alpha)/(2*(2-alpha^2)));

    apb = ab.*(d - ab + alpha*bb );
    bpb = bb.*(d - bb + alpha*ab );

    figure
    subplot(2,2,1)
    hold on
    plot(d,aa,'k')
    plot(d,ab,'r')
    hold off
    set(gcf,'Color',[1 1 1])
    title('Firm A Price')
    xlabel('d')
    ylabel('Price')
    
    subplot(2,2,2)
    hold on
    plot(d,ba,'k')
    plot(d,bb,'r')
    hold off
    set(gcf,'Color',[1 1 1])
    title('Firm B Price')
    xlabel('d')
    ylabel('Price')
    
    subplot(2,2,3)
    hold on
    plot(d,apa,'k')
    plot(d,apb,'r')
    hold off
    set(gcf,'Color',[1 1 1])
    title('Firm A Profits')
    xlabel('d')
    ylabel('Profits')
    
    subplot(2,2,4)
    hold on
    plot(d,bpa,'k')
    plot(d,bpb,'r')
    hold off
    set(gcf,'Color',[1 1 1])
    title('Firm B Profits')
    xlabel('d')
    ylabel('Profits')
    legend('Part A','Part B','Location','NorthWest')
    
    saveas(gcf,['question6' num2str(count) '.png'])
end