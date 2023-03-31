%% ME 674 - Homework 7
clear, clc, close all

%% Setup

%defining grid variables 
a = 5;  %Defines how large your view window is
n = 500; % number of intervals (Smoothes the lines)

X = linspace(a, -a, n);
Y = linspace(a, -a, n);

[x,y]=meshgrid(X,Y);

%Defining complex variable z
z = x + i*y;

%   Take note: top-level code is in cartesian coords. Functions describing
%       a flow will switch to polar if the equation warrants it.

%% Plotting
%Define your desired compound flow here by strining together multiple flow
%   types. Check below in the functions section to see what is available

%   Lifting Cylinder
F1 = uniformFlow(1, 1, z) + doublet(1.5, z, 0, 0)  + irrVortex(-3, z, 0, 0);

figure (1)
contour(x,y,imag(F1),70,'Linewidth',1)
axis square
title("Lifting Cylinder")
xlabel("X Position")
ylabel("Y Position")

% %   Source in Uniform Flow
% F2 = sourceSink(10, z, 0, 0) + uniformFlow(1, 1, z);
% 
% figure (2)
% contour(x,y,imag(F2),22,'Linewidth',1)
% axis square
% title("Source in Uniform Flow")
% xlabel("X Position")
% ylabel("Y Position")

% %   Two Sources
% F3 = sourceSink(10, z, -10, 0) + sourceSink(10, z, 10, 0);
% 
% figure (3)
% contour(x,y,imag(F3),25,'Linewidth',1)
% axis square
% title("Two Sources Interacting")
% xlabel("X Position")
% ylabel("Y Position")

% %   Attempt at a NACA 2412 Airfoil
% data = readmatrix('NACA_2412.txt').*20;
% x_data = data(:,1)-3;
% y_data = data(:,2);
% F4 = uniformFlow(2, 1, z);
% 
% for i = 1:length(x_data)
%     F4 = F4 + irrVortex(-1.5, z, x_data(i), y_data(i)) + doublet(0.3, z, x_data(i), y_data(i));
% end
% 
% figure (4)
% contour(x,y,imag(F4),45,'Linewidth',1)
% axis square
% title("Attempt at NACA 2412 Airfoil via Vortex Sheet")
% xlabel("X Position")
% ylabel("Y Position")

% %   Array of Lifting Cylinders
% s = 7;
% v = 20;
% F5 = doublet(s, z, -5, 0)  + irrVortex(-v, z, -5, 0) + uniformFlow(3, 1, z)...
%     + doublet(s, z,-2.2, 3)  + irrVortex(-v, z, -2.2, 3)...
%     + doublet(s, z, 2.2, 3)  + irrVortex(-v, z, 2.2, 3)...
%     + doublet(s, z, -2.2, -3)  + irrVortex(-v, z, -2.2, -3)...
%     + doublet(s, z, 2.2, -3)  + irrVortex(-v, z, 2.2, 3)...
%     + doublet(s, z, 5, 0)  + irrVortex(-v, z, 5, 0);
% 
% figure (5)
% contour(x,y,imag(F5),25,'Linewidth',1)
% axis square
% title("Array of Lifting Cylinders")
% xlabel("X Position")
% ylabel("Y Position")
% 
% %   Imbalanced Source/Sink Pair
% F6 = sourceSink(5, z, -10, 0) + sourceSink(-10, z, 10, 0);
% 
% figure (6)
% contour(x,y,imag(F6),25,'Linewidth',1)
% axis square
% title("Imbalanced Source/Sink Pair")
% xlabel("X Position")
% ylabel("Y Position")

%#######    Neither of these flows are working at the moment    #########

% %   Weird Flows
% F5 = sourceSink(1, z, 0, 0) + irrVortex(-1, z, 0, 0);
% 
% figure (5)
% contour(x,y,imag(F5),25,'Linewidth',1)
% axis square
% title("Source with Vortex")
% xlabel("X Position")
% ylabel("Y Position")

% %   Weird Flows
% F5 = sourceSink(1, z, -10, 5) + sourceSink(1, z, 10, 5) + sourceSink(-1, z, -5, 0);
% 
% figure (5)
% contour(x,y,imag(F5),25,'Linewidth',1)
% axis square
% title("Source with Vortex")
% xlabel("X Position")
% ylabel("Y Position")


%% Functions
function func_uniformFlow = uniformFlow(U, n, z)
%Plots an irrotational vortex using complex potential
    %U is velocity of uniform flow
    %n is used to determine the flow angle from the equation
    %   alpha = pi/n (n = 2 is flow impinging on a wall)
    %z carries the meshgrid data
    
    func_uniformFlow = U.*(z^n);

end

function func_ss = sourceSink(q, z, x, y)
%Plots a source or sink using complex potential
    %x, y are center of the source or sink
    %q is the strength of the source or sink
    %z carries the meshgrid data

%     z_prime = x + y.*i;
%     
    z(abs(z)<1) = 0.001;
%     
%     func_ss = (q/(2*pi)).*log(z - z_prime);

    %Defining polar variables
    diff_x = real(z) - x;
    diff_y = imag(z) - y;
    
    r = sqrt(diff_x.^2 + diff_y.^2);
    theta = atan2(diff_y,diff_x);
    
%     contourf(theta)
    
    z = r.*exp(i.*theta);

    func_ss = (q/(2*pi)).*log(r) + i.*theta.*(q/(2*pi));
    
%Junk pile (keep for posterity and record-keeping)
% theta = mod(atan2(diff_y,diff_x), 2*pi);
    % increasees atan2's range to 0 to 2pi, but does not change our
    % discontinuity error. Flips the weird graph effects to 0 as opposed to
    % pi, so that's nice
    
%theta = unwrap(mod(atan2(diff_y,diff_x), 2*pi));
    %Changes nothing except directs the weird behavior downward
    


end

function func_ir = irrVortex(G, z, x, y)
%Plots an irrotational vortex using complex potential
    %x, y are center of vortex
    %G (Gamma) is the strength of circulation of vortex
    %z carries the meshgrid data
    
%     z_prime = x + y*i;
%     
%     func_ir = ((-i*G)/(2*pi)).*log(z - z_prime);

    diff_x = real(z) - x;
    diff_y = imag(z) - y;
    
    r = sqrt(diff_x.^2 + diff_y.^2);
    theta = atan2(diff_y,diff_x);
    
    z = r.*exp(i.*theta);

    func_ir = (G/(2*pi)).*theta - i.*log(r).*(G/(2*pi));

end

function func_doublet = doublet(d, z, x, y)
%Plots an irrotational vortex using complex potential
    %x, y are center of the doublet
    %d is the strength of the doublet
    %z carries the meshgrid data
    
    z_prime = x + y*i;
    
    func_doublet = d./((2*pi).*(z-z_prime));

end
