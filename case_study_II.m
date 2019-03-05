%{
    Title: Heat Transfer Case Study II HDE Solver
    Author: James Emerson Parkus
    Date: 3/5/19
    Purpose: Solve the HDE for a 1 dimensional solid with adiabatic and
    convective boundary conditions. The solution output is a temperature
    profile and a temperature at the centerline of the solid.
%}

clc
clear
close all

load coefficients

%% Constant
k = 0.04; % [W/(m K)] - Thermal conductivity of perlite (ranges from 0.04 -> 0.06)
cp = 837; % [J/(kg K)] - Specific heat for constant pressure
rho = 2200; % [kg/m^3] - Density (ranges from 2200 -> 2400)
h = 30; % [W/(m^2 K)] - Convective coefficient for air during fire
L = 0.0508; % [m] - Length of safe wall

T_inf = 1088.706; % [K] - Temperature of the outside air temperature during fire
T_i = 293.15; % [K] - The initial temperature of the safe before the fire

%% Analysis
alp = k/(cp*rho); % [m^2/s] - Thermal diffusivity
L_c = L/2; % [m] - Characteristic length
Bi = h*L_c/k; % [-] - Biot number
lam = 1.3; % [-] - Solution coefficient
C_1 = 4*sin(lam)/(2*lam + sin(2*lam)); % [-] - Coefficient C_1 for solution

i = 1;
j = 1;
for t = 3000:100:6000
    Fo = alp*t/L_c^2;
    
    for x = 0:0.001:1
%         temp = F(x,Fo,lam,C_1);
        T = T_inf + (T_i - T_inf)*C_1*exp(-lam^2*Fo)*cos(lam*x);
        
        temperature(j,i) = T;
        distance(j,i) = x;
        j = j + 1;
    end
    
    fourier_number(i) = Fo;
    time(i) = t;
    i = i + 1;
    j = 1;
end

hold on
grid on
for i = 1:1:1001

plot(distance(:,i),temperature(:,i));

end

%% Functions
% 
% function T = F(x,Fo,lam,C_1)
% global T_inf T_i
% 
% T = T_inf + (T_i - T_inf)*C_1*exp(-lam^2*Fo)*cos(lam*x);
% end









