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

%% Given
k_perlite = 0.06; % [W/(m K)] - Thermal conductivity of perlite (ranges from 0.04 -> 0.06)
k_plastic = 0.27; % [W/(m K)] - Thermal conductivity of plastic wall
cp = 837; % [J/(kg K)] - Specific heat for constant pressure
rho = 2200; % [kg/m^3] - Density (ranges from 2200 -> 2400)
h = 30; % [W/(m^2 K)] - Convective coefficient for air during fire
L_perlite = 0.0239; % [m] - Length of safe wall
L_plastic = 0.00127; % [m] - Length of plastic wall (assumed to be nylon 6/6)

T_inf = 1116.483; % [K] - Temperature of the outside air temperature during fire
T_i = 293.15; % [K] - The initial temperature of the safe before the fire

%% Analysis
R_eff_prime = 1/h + L_plastic/k_plastic; % [K/W] - Effective thermal resistance of convection and conduction on plastic wall
U = 1/R_eff_prime; % [W/K] - Transforming R_eff_prime into a form that can be used for the HDE solution
alp = k_perlite/(cp*rho); % [m^2/s] - Thermal diffusivity
L_c = L_perlite; % [m] - Characteristic length
Bi = U*L_c/k_perlite; % [-] - Biot number

for i = 2:1:length(coefficients)
    if Bi > coefficients(i-1,1) && Bi < coefficients(i,1)
        if abs(Bi - coefficients(i-1,1)) > abs(Bi - coefficients(i,1))
            lam = coefficients(i,2); % [-] - Solution coefficient
        else
            lam = coefficients(i-1,2);
        end
    end
end

C_1 = 4*sin(lam)/(2*lam + sin(2*lam)); % [-] - Coefficient C_1 for solution

i = 1;
j = 1;
for t = 3000:100:6000
    Fo = alp*t/L_c^2;
    
    for x = 0:0.001:1
        T = T_inf + (T_i - T_inf)*C_1*exp(-lam^2*Fo)*cos(lam*x);
        
        temperature(j,i) = T;
        distance(j,i) = x;
        j = j + 1;
    end
    
    fourier_number(i,1) = Fo;
    time(i,1) = t;
    i = i + 1;
    j = 1;
end

hold on
grid on
for i = 1:1:1
    plot(distance(:,i),temperature(:,i));
    title('HDE Solution');
    xlabel('x^* [-]');
    ylabel('Temperature [K]');
end
