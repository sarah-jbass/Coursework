clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 711\Lones Smith\Problem Sets\PS1'

%Written by Sarah Bass, 11/2/20
grid = 0:0.0001:1;

BR1 = grid<=0.6;
BR2 = grid<=(8/11);
figure
hold on
plot(BR1,grid,'r-')
plot(grid,BR2,'b-')
plot([0 1 8/11],[1 0 0.6],'ko')
hold off
set(gcf,'Color',[1 1 1])
xlim([0 1.1])
ylim([0 1.1])
ylabel('\sigma_2')
xlabel('\sigma_1')
legend('Best response 1','Best response 2','Equilibria')
title('Best response curves')
saveas(gcf,'best_response.png')