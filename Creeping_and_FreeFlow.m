%ME 674 Homework 5 - Problem 4
clc

%% Generating PLot Data
%All values that can be are set to unity unless otherwise noted
U = 1;      %Freestream velocity
a = 1;

b = a*2.5;
c = -a*2.5;
n = 60; % number of intervals (Smoothes the lines)

[x,y]=meshgrid((c:(b-c)/n:b),(c:(b-c)/n:b)');

%Need to remove data computed when r<a and set it to zero
for i = 1:length(x)
   for k = 1:length(x)
      if sqrt(x(i,k).^2 + y(i,k).^2) < a
         x(i,k) = 0;
         y(i,k) = 0;
      end
   end
end

%Coordinate switch to polar
r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);

% Freestream profile
psi_free = (1/2).*U.*(r.^2).*(1-((a.^3)./(r.^3))).*(sin(theta)).^2; %pg 326, K&C, eq 7.85

%Creeping Profile
psi_creep = (U.*(r.^2).*(sin(theta)).^2).*((1/2)-((3.*a)./(4.*r))+((a.^3)./(4.*r.^3)));  %pg 444, K&C, eq 9.53
%psi_creep = (U.*(r.^2).*(sin(theta)).^2)

%% Plotting stream functions
% Freestream line Plot
figure(1)
contour(x,y,psi_free,32,'Linewidth',1)
hold on 
axis square
title('\psi_{free}')
circle(0, 0, a);
hold off

% Creeping line Plot
figure(2)
contour(x,y,psi_creep,32,'Linewidth',1)
hold on 
axis square
title('\psi_{creeping}')
circle(0, 0, a);
hold off


%% Functions 
function drawCircle = circle(x,y,r)
hold on
theta = 0:pi/50:2*pi;
xunit = r*cos(theta) + x;
yunit = r*sin(theta) + y;
drawCircle = plot(xunit, yunit, 'Linewidth',3);
hold off
end