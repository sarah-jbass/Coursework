% Econ 711 Problem Set 1
%
% Prepared by Sarah Bass
% Last updated: 09/10/15
pause('on')

clear;
clc;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 711\Problem Sets\PS1'

y = [-20 40; -40 70; -70 90; 0 0];
prices = [7 4; 5 5; 4 8];
[line,itsct] = do_lines(y(1:end-1,1),y(1:end-1,2),prices);

%grid = -100:0;
grid = -100:25;
line1 = line(1,1)*grid + line(1,2);
line2 = line(2,1)*grid + line(2,2);
line3 = line(3,1)*grid + line(3,2);
vecscale = 1;

figure
hold on
h = fill([grid(end) itsct(:,1)' grid(1) grid(1) grid(end)],[ line1(end) itsct(:,2)' line3(1) -25 -25],[0.95 0.95 0.95]);
h(1).EdgeColor = [0.95 0.95 0.95];
plot(grid,line1,'k')
plot(grid,line2,'k')
plot(grid,line3,'k')
for i=1:3
plot([y(i,1) y(i,1)+vecscale*prices(i,1)],[y(i,2) y(i,2)+vecscale*prices(i,2)],'b-')
end
%plot(grid,yg,'r');
plot(y(1:end-1,1),y(1:end-1,2),'r*')
xline(0)
yline(0)
hold off
set(gcf,'Color',[1 1 1])
title('Largest production set that can rationalize the data')
%legend('production function','Location','NorthEast')
xlim([-100 25])
ylim([-25 100])
saveas(gcf,'largestprod.png') %savefig


figure
hold on
h1 = fill([0 y(1:end-1,1)' -100 -100 0],[0 y(1:end-1,2)' y(end-1,2) -25 -25],[0.95 0.95 0.95]);
h1(1).EdgeColor = [0.95 0.95 0.95];
pause(3)
plot([0; y(1:end-1,1)],[0; y(1:end-1,2)],'k') % Need to make this convex!
pause(3)
plot(y(1:end,1),y(1:end,2),'r*')
xline(0)
yline(0)
hold off
set(gcf,'Color',[1 1 1])
title('Free disposal, convex, shutdown property')
%legend('production function','Location','NorthEast')
xlim([-100 25])
ylim([-25 100])
saveas(gcf,'freedispshutdown.png')


function [line,itsct] = do_lines(x,y,prices)
% Calculates lines perp to prices going through the pts (takes first col of 
%vector y as x and second column of vector y as y)
line = [0*x 0*y];
itsct = zeros(size(line,1)-1,size(line,2));
for i=1:length(x)
    line(i,1) = -(prices(i,1)/prices(i,2)); %first col is the slope
    line(i,2) = y(i) - line(i,1)*x(i); % second col is the intercept
end
for i=1:length(x)-1
    itsct(i,1) = (line(i+1,2) - line(i,2))/(line(i,1) - line(i+1,1)); %x coordinate of intersection
    itsct(i,2) = line(i,1)*itsct(i,1) + line(i,2); %y coordinate of intersection
end
end