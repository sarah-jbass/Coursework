clear;
clc;
close all;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Spring 2021\Econ 714\Dmitry Mukhin\Problem Sets\PS1'

%Written by Sarah Bass, 1/26/21

%Params
kgrid = 150:0.5:200; % capital grid for phase diagram
cgrid = 0:0.0001:5;
D0 = 0;
D1 = 1;
T = 12;
ns = 120;
psigma = 1;
palpha = 1/3;
pbeta = 0.99^(1/12);
pdelta = 0.01;
F = @(K) K.^palpha; %why is there no A?
Fp = @(K) palpha*K.^(palpha - 1);
U = @(C) ( C.^(1-psigma)- 1)./(1-psigma);
Up = @(C) (C.^(-psigma));

%Calculate Steady State
deltak = F(kgrid) - pdelta*kgrid - D0; 
% delta c line intersects kss
kss = (((1/pbeta) - 1 + pdelta)/palpha)^(1/(palpha - 1));
css = F(kss) - pdelta*kss - D0;
[kcheck,ccheck] = iterate(kss,css, D0, psigma, palpha, pbeta, pdelta);
if abs(kcheck - kss)+abs(ccheck-css)>1e-10; warning('failed check, ss is not stationary'); end

%Q2: Calculate Original Saddle Path
recalc=0;
if recalc==1
    saddlepath = saddle(kgrid,cgrid,kss,css,D0,psigma,palpha,pbeta,pdelta,ns); %(kgrid, cgrid, kss, css, D, psigma, palpha, pbeta, pdelta,ns)
    save 'saddle_data' 'saddlepath'
else
    load 'saddle_data'
end

%Q3: Calculate Transition Path
ns = 600;
D = D0 + zeros(1,ns);
D(T+1) = D1;
recalc2=0;
if recalc2==1
    [Ktraj, Ctraj, diff] = trajectory(kss,css,D,palpha,psigma,pbeta,pdelta);
    save 'traj_data' 'Ktraj' 'Ctraj'
else
    load 'traj_data'
end

%Q2: Phase Diagram (Zoomed In)
figure
hold on 
plot(kgrid, deltak, 'k')
plot([kss kss],[0 5], 'b')
plot(kgrid, saddlepath, 'r')
hold off
set(gcf,'Color',[1 1 1])
title('Phase Diagram')
xlabel('K')
ylabel('C')
ylim([3,5])
legend('\Delta K = 0', '\Delta C = 0', 'Saddle Path', 'Location', 'SouthEast')
annotation('arrow',[0.3 0.3],[0.2 0.25])
annotation('arrow',[0.3 0.35],[0.2 0.2])
annotation('arrow',[0.3 0.3],[0.8 0.85])
annotation('arrow',[0.3 0.25],[0.8 0.8])
annotation('arrow',[0.75 0.75],[0.35 0.3])
annotation('arrow',[0.75 0.8],[0.35 0.35])
annotation('arrow',[0.7 0.7],[0.8 0.75])
annotation('arrow',[0.7 0.65],[0.8 0.8])
saveas(gcf,'q2_phasediagram.png')

%Q2: Phase Diagram (Zoomed Out)
kgrid2 = 0:1200;
delk2 = F(kgrid2) - pdelta*kgrid2 - D0; 
ii = delk2>=0;
kgrid2 = kgrid2(ii);
delk2 = delk2(ii);

figure
hold on 
plot(kgrid2, delk2, 'k')
plot([kss kss],[0 5], 'b')
plot(kss,css,'r*')
plot(0,0,'r*')
plot(kgrid2(end),delk2(end),'r*')
hold off
set(gcf,'Color',[1 1 1])
title('Phase Diagram')
xlabel('K')
ylabel('C')
ylim([0,5])
legend('\Delta K = 0', '\Delta C = 0', 'Equilibria', 'Location', 'NorthEast')
saveas(gcf,'q2_phasediagram2.png')

t=(-99:ns);
%Q3: Capital Trajectory
figure
hold on
plot(t, [repmat(kss,1,100) Ktraj],'k')
plot(t,kss+ t*0,'b-')
plot([0 0],[168 172],'Color',[0 0.6 0])
plot([12 12],[168 172],'r-')
hold off
set(gcf,'Color',[1 1 1])
title('Capital Trajectory Over Time')
xlabel('time')
ylabel('K')
legend('Capital Trajectory', 'Steady State Capital', 'Time of News', 'Time of Shock','Location', 'SouthEast')
saveas(gcf,'q3_ktraj.png')

%Q3: Consumption Trajectory
figure
hold on
plot(t, [repmat(css,1,100) Ctraj],'k')
plot(t,css+ t*0,'b-')
plot([0 0],[3.82 3.85],'Color',[0 0.6 0])
plot([12 12],[3.82 3.85],'r-')
hold off
set(gcf,'Color',[1 1 1])
title('Consumption Trajectory Over Time')
xlabel('time')
ylabel('K')
ylim([3.82, 3.85])
legend('Consumption Trajectory', 'Steady State Consumption', 'Time of News', 'Time of Shock','Location', 'SouthEast')
saveas(gcf,'q3_ctraj.png')

%Q3: Phase Diagram
figure
hold on
plot(kgrid, deltak, 'k')
plot([kss kss],[0 5], 'b')
plot([kss Ktraj],[css Ctraj],'r')
plot(Ktraj(1),Ctraj(1),'mo')
plot(Ktraj(T),Ctraj(T),'mx')
plot(Ktraj(T+1),Ctraj(T+1),'m+')
plot(kss,css,'m*')
hold off 
set(gcf,'Color',[1 1 1])
title('Phase Diagram')
xlabel('K')
ylabel('C')
ylim([3.78, 3.9])
xlim([169.6, 170.8])
legend('\Delta K = 0', '\Delta C = 0', 'Trajectory', 'Initial Jump', 'Position before shock', 'Position at shock', 'Final position (steady state)', 'Location', 'NorthWest')
saveas(gcf,'q3_phasediagram.png')