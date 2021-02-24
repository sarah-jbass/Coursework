clear;
clc;

cd 'C:\Users\19195\OneDrive\Documents\University of Wisconsin\First Year\Fall Classes\Econ 712\Problem Sets\PS3'

% Written by Sarah Bass

%%%%%%%%%%%%Part A%%%%%%%%%%%%

%Creating a matrix with p1/p2, and corresponding c1, c2, and U
A = zeros(4,300);
%p1/p2:
A(1,:) = 0.05:0.05:15;
%c1:
for i = 1:300
    A(2,i) = 10/(8 + 2*(A(1,i)^2));
end
%c2:
for i = 1:300
    A(3,i) = 2 - A(1,i)*A(2,i);
end
%U_p:
for i = 1:300
    A(4,i) = 10*A(2,i) - 4*A(2,i)^2 + 4*A(3,i)- A(3,i)^2;
end

c_1 = -sqrt(41/16) + (5/4):0.01:sqrt(41/16) + (5/4);
c_2 = zeros(6,321);
a=1;
for i = [10.1538461538462 9.88235294117647 9]
    c_2(a,:) = sqrt((41/4) - i - 4*(c_1-(5/4)).^2) +2; %top of ellipse
    a=a+1;
    c_2(a,:) = -sqrt((41/4) - i - 4*(c_1-(5/4)).^2) +2; %bottom of ellipse
    a=a+1;
end

%Creating lines for our budget constraints
bc_1= -0.25*c_1 + 2;
bc_2= -0.5*c_1 + 2;
bc_3= -1*c_1 + 2;

%Plotting figure
figure
hold on;
for i = 1:6
   plot (c_1,c_2(i,:),'b');
end
plot (A(2,:), A(3,:), 'r');
plot (c_1,bc_1,'k');
plot (c_1,bc_2,'k');
plot (c_1,bc_3,'k');
plot (5/4,2,'r*');
axis([-0.25 2.5 0 4])
xlabel('c_1')
ylabel('c_2')
title('Utility, budget constraints, and the offer curve')
saveas(gcf,'Question2_A.png')

%%%%%%%%%%%%Part B%%%%%%%%%%%%

%Creating line segments
c1_0 = [0 0];
c2_0 = [2 3];
c1_a = [0 2/3];
c2_a = [2 2/3];
c1_b = [2/3 1/3];
c2_b = [2/3 1/3];
c1_c = [1/3 1];
c2_c = [1/3 0];

%Plotting figure
figure
hold on;
plot (c1_0,c2_0,'r');
plot (c1_a,c2_a,'r');
plot (c1_b,c2_b,'r');
plot (c1_c,c2_c,'r');
plot (1/3,1/3,'r*');
plot (2/3,2/3,'r*');
axis([0 2 0 3])
xlabel('c_1')
ylabel('c_2')
title('Offer Curve')
saveas(gcf,'Question2_B.png')

%%%%%%%%%%%%Part B%%%%%%%%%%%%

%Creating line segments
c1_a = [0 4];
c2_a = [12 4];
c1_b = [4 7];
c2_b = [4 7];
c1_c = [7 21];
c2_c = [7 0];

%Plotting figure
figure
hold on;
plot (c1_a,c2_a,'r');
plot (c1_b,c2_b,'r');
plot (c1_c,c2_c,'r');
plot (4,4,'r*');
plot (7,7,'r*');
axis([0 25 0 15])
xlabel('c_1')
ylabel('c_2')
title('Offer Curve')
saveas(gcf,'Question2_C.png')