% pause('on')

clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Problem Sets\PS7'

%Written by Sarah Bass, 10/22/20

%Parameters
y=1;
pbeta = 0.999;
palpha = 0.99999;
pt = 0:0.01:y-0.01;
pt1 =pt./(pbeta*(y - pt)) - palpha/pbeta;

%Graph
figure 
hold on;
plot (pt, pt1, 'r')
plot([-100 100],[0 0],'k--')
xlabel('P_t')
ylabel('P_{t+1}')
xlim([0 1.2])
ylim([-5 20])
set(gca,'xtick',[],'ytick',[])
title('Law of Motion for House Prices in Equilibrium')
saveas(gcf,'law_of_motion.png')