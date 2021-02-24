
clear;
clc;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 711\Problem Sets\PS2'

% Written by Sarah Bass

y_1 = -100:0.01:0; %the 0.01 is the step between elements
y_2 = (-1*y_1).^(2/3); %. means that it's element wise

figure
hold on;
plot (y_1,y_2,'b')
plot([0 0],[-100 100],'k')
axis([-100 20 -5 25])
p = patch([y_1 0 -100],[y_2 -100 -100], 'b');
p.FaceColor = [0.9 0.9 0.9];
plot([-100 100],[0 0],'k')
xlabel('y_1')
ylabel('y_2')
title('Production set with B=1')
saveas(gcf,'Question_2a.png')
