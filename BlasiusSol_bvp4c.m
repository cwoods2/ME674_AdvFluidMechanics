%% Problem 2 - Numerical Solution to the Blasius Equation (bvp4c)
% f''' + (1/2)f*f'' = 0
%   where f = f(eta)
%   f'(0) = 0 and f'(inf) = 1
clear, clc, close all

%Initializing bvp
eta = linspace(0,30,500);
init = bvpinit(eta, @guess);

%Solving using bvp4c
sol = bvp4c(@Blasius_bvp4c, @Blasius_bc, init);

%Extracting data of interest
f = transpose([sol(:,1).y]);
eta = transpose([sol(:,1).x]);

table(:,1) = eta;
table(:,2) = f(:,1);
table(:,3) = f(:,2);
table(:,4) = f(:,3);

%Plotting results
figure(1)
hold on
plot(eta, f, "Linewidth", 2)
yline(1,'--')
axis([0 7 0 3])
title("Solutions to the Blasius Equation")
xlabel("\eta")
ylabel("Function Values")
legend({"f", "{\partialf}/_{\partial\eta}", "{\partial^2f}/_{\partial\eta^2}"}, 'Location','northwest')
hold off

figure(2)
hold on
plot(eta, f(:,2), "Linewidth", 2)
yline(1,'--')
axis([0 7 0 1.1])
title("Blasius Similarity Solution (Figure 10.5 from K&C)")
xlabel("\eta")
ylabel("u/U")
hold off


%% Problem 3 - Redoing Example 10.3
f1 = f(:,1);
psi = [0.001, 0.004, 0.01, 0.02];
U = 1;                  %Freestream velocity [m/s]
nu = 1.46e-5;           %Kinematic Viscosity of air [m^2/s]

%Calculating values for the given constant streamline values
psi_x = (psi.^2)./((f1.^2)*nu*U);
psi_y = (eta.*psi)./(f1)*U;

%Defining delta 99
x = linspace(0,3, 200);
delta_99 = 4.92*sqrt((nu.*x)./U);

%Plotting results
figure(3)
hold on
plot(psi_x, psi_y)
plot(x, delta_99, 'r', 'Linewidth', 2)
axis([0 3 0 0.04])
title("Selected Constant Freestram Profiles with \delta_{99} included")
xlabel("x [m]")
ylabel("y [m]")
legend("\psi = 0.001", "\psi = 0.004", "\psi = 0.01", "\psi = 0.02", "\delta_{99}", 'Location','northwest')
hold off

%% Functions Used for HW6 Part 2 - Blasius Solution
function func = Blasius_bvp4c(eta,f) %Blasius Function 

    func = [f(2); f(3); -0.5*f(1)*f(3)];
end

function res = Blasius_bc(y1,y2) % boundary conditions
    res = [y1(1); y1(2); y2(2)-1];
end

function g = guess(x)   %Initial Guess
    g = [0, 0, 0];
end
