%% ME 674 - Problem Set 8
%Camden Woods, Spring 2022
clear, clc, close all

%% Setup
%Setting up meshgrid environment
q = 4;  %Defines how large your view window is
n = 500; % number of intervals (Smoothes the lines)

%zed plane setup
X = linspace(q, -q, n);
Y = linspace(q, -q, n);

[x,y]=meshgrid(X,Y);

%Defining complex variable z
zeta = x + i.*y;

%Changing contour bounds
contour_amounts = [-10:0.5:10];

%% General Problem 1
%Plotting a plane wall using method of images
wall_Flow = planeWall(10, zeta, 0, 1);

% figure (1)
% hold on
% contour(x,y,imag(wall_Flow),contour_amounts,'Linewidth',1)
% axis square
% title("Horizontal Plane Wall From Two Counter-Rotating Vortices")
% xlabel("x Position")
% ylabel("iy Position")
% hold off

%% General Problem 2
%Generating Different Airfoils and Plotting Lift Characteristics
c = 1;

%Airfoil Coords
xpos = -0.17; 
ypos = 0.15;

% xpos = -0.1; %High camber airfoil
% ypos = 0.15;    
% 
% xpos = -1.5; %Flat plate coords
% ypos = 0;

%Calculating the x and y offsets
%   With this implementation, you either get to specify the radius of the circle
%   or the offsets. Defining both will cause the circle to fall off of point c
r = (ypos^2 + (xpos - c)^2)^0.5;  %For all other shapes
beta_prime = asind(ypos./r);
% r = 2;    %For ellipses

%Eliminating interior flow
zeta_0 = xpos + i*ypos; %Accounts for moving the zeta-circle from the default position (the origin)

r_zeta = abs(zeta-zeta_0);
zeta(r_zeta <= r) = nan;

%Defining zeta space terms
%           cylinderFlow(vel, a, alpha, zeta, x, y, Beta)
vel = 1;
alpha = 20;

zeta_flow = cylinderFlow(vel, r, alpha, zeta, xpos, ypos, beta_prime);
zeta_circle = circle(r, xpos, ypos);

%Transforming terms to zed-space
airfoil = zeta_circle + (c^2./zeta_circle);

r_zeta = abs(zeta);
zeta(r_zeta <= r) = nan;

z = zeta + ((c^2)./zeta);

%Computing Lift values
rho = 1.225;    %Density of air

    %Calculating components of velocity
    U = vel*cosd(alpha);
    V = vel*sind(alpha);
    
    U_inf = sqrt(U^2 + V^2);
    
    %Calculating Gamma (for cambered airfoils)
    Gamma = 4*pi*U_inf*r*sind(alpha + beta_prime);
    
Lift = U_inf*rho*Gamma;

%Computing Pressure Values
P_SL = 101325;  %Pa
    
    %Complex Velocity, W, and velocity in z direction
    W = (U - i*V) - (U + i*V).*((r.^2)./(zeta-zeta_0).^2) + (i*Gamma.*r)./((2*pi).*(zeta-zeta_0));
    
    df_dz = W.*((zeta.^2)./(zeta.^2 - 1));
    
    %Calculating Pressure from Bernoulli
    P = P_SL + 0.5*rho*(df_dz).^2;
    
    %Calculating Pressure Coefficient
    Cp = (P - P_SL)./(0.5*rho*U_inf^2);

figure (2)
subplot(1,2,2)
hold on
contour(real(z),imag(z),imag(zeta_flow),contour_amounts,'Linewidth',1)
plot(real(airfoil), imag(airfoil), 'k', 'Linewidth', 2);
xline(0, '--');
yline(0, '--');
axis equal
title("Transformed Flow in Zed-Plane")
subtitle("Lift per unit span: " + Lift + "  [N/m], Angle of attack: " + alpha + " [deg]")
xlabel("X Position")
ylabel("iY Position")
hold off
subplot(1,2,1)
hold on
contour(x,y,imag(zeta_flow),contour_amounts,'Linewidth',1)
plot(real(zeta_circle), imag(zeta_circle), 'b', 'Linewidth', 2);
circle_reference(c, 0, 0);
xline(0, '--');
yline(0, '--');
axis square
title("Cylinder in Angled Flow Using the Circle Theorem, Zeta-Plane")
xlabel("\xi Position")
ylabel("i\eta Position")
hold off

figure (3)
subplot(1,2,1)
hold on
contourf(real(z), imag(z), abs(P))
plot(real(airfoil), imag(airfoil), 'k', 'Linewidth', 2);
axis equal
title("Pressure Distribution on Airfoil")
xlabel("X Position")
ylabel("Y Position")
hold off
subplot(1,2,2)
hold on
contourf(real(z), imag(z), abs(Cp))
plot(real(airfoil), imag(airfoil), 'k', 'Linewidth', 2);
axis equal
title("Coefficient of Pressure Distribution on Airfoil")
xlabel("X Position")
ylabel("Y Position")
hold off

%% Functions
function func_wall = planeWall(G, z, x, y)
%Plots a plane wall using method of images
%coordinates x, y, denote the location of the center of the wall

    %Standard Irrotational Vortex
    z_prime = x + y*i;
    
    func_irr = ((-i*G)/(2*pi)).*log(z - z_prime);
    
    %Conjugate Irrotational Vortex
    z_prime_conj = x - y*i;
    
    func_irr_conj = ((i*G)/(2*pi)).*log(z - z_prime_conj);

    func_wall = func_irr + func_irr_conj;

end

function func_flowCyl = cylinderFlow(vel, a, alpha, zeta, x, y, beta)
%Plots a cylinder in a uniform flow angled at a given angle of attack
%a is the radius of the cylinder
%alpha is the angle of attack of the flow, input in degrees
%U is the freestream strength of the flow

    %Calculating components of velocity
    U = vel*cosd(alpha);
    V = vel*sind(alpha);
    
    U_inf = sqrt(U^2 + V^2);
    
    %Implementing transformation
    zeta_prime = x + i*y;
    
    %Calculating Gamma (for cambered airfoils)
    Gamma = 4*pi*vel*a*sind(alpha+beta);

    func_flowCyl = ((U - i*V).*(zeta-zeta_prime)) + (U + i*V).*((a^2)./(zeta-zeta_prime)) + (i*Gamma/(2*pi))*log((zeta-zeta_prime)./a);
end

function h = circle_reference(r, x, y)
    hold on
    theta = 0:pi/50:2*pi;
    xunit = r*cos(theta) + x;
    yunit = r*sin(theta) + y;
    h = plot(xunit, yunit, '--', 'Linewidth', 2);
    hold off
end

function func_Circle = circle(r, x, y)
    theta = 0:pi/50:2*pi;
    xunit = r*cos(theta) + x;
    yunit = r*sin(theta) + y;
    func_Circle = xunit + i.*yunit;
end