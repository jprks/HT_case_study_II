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

%% Given
global cp rho
k_perlite = 0.0775; % [W/(m K)] - Thermal conductivity of perlite (ranges from 0.04 -> 0.06)
k_plastic = 0.46; % [W/(m K)] - Thermal conductivity of plastic wall
cp = 837; % [J/(kg K)] - Specific heat for constant pressure
rho = 1200; % [kg/m^3] - Density (ranges from 2200 -> 2400)
h = 30; % [W/(m^2 K)] - Convective coefficient for air during fire
L_perlite = 0.0239; % [m] - Length of safe wall
L_plastic = 0.00127; % [m] - Length of plastic wall (assumed to be nylon 6/6)

T_inf = [1116.483;404.5]; % [K] - Temperature of the outside air temperature during fire
T_s_i = [293.15;294.45]; % [K] - The initial temperature of the safe before the fire

%% Uncertainties
del_k_Perl = 0.0125; % [W/(m K)] - Uncertainty of Perlite k, +/- 0.0125
del_k_HDPE = 0.04; % [W/(m K)] - Uncertainty of HDPE k, +/- 4.00
del_h = 3.00; % [W/(m^2 K)] - Uncertainty of Air h, +/- 3.00
del_L = 0.0005*0.0254; % [m] - Uncertainty of Length l, +/- 0.0005"
% del_rho = 100; % [kg/m^3] - Uncertainty of Density rho
% +/- 5% for Constant Specific Heat

for n = 1:1:10
    pungent = n/10;
    dev_k_Perl(n,1) = del_k_Perl*pungent;
    dev_k_HDPE(n,1) = del_k_HDPE*pungent;
    dev_h(n,1) = del_h*pungent;
    dev_L(n,1) = del_L*pungent;
end

for i = 1:1:10
    k_Perl(i,1) = k_perlite - dev_k_Perl(11-i,1);
    k_HDPE(i,1) = k_plastic - dev_k_HDPE(11-i,1);
    h_Air(i,1)  = h - dev_h(11-i,1);
    L_Perl(i,1) = L_perlite - dev_L(11-i,1);
    L_HDPE(i,1) = L_plastic - dev_L(11-i,1);

    k_Perl(i+11,1) = k_perlite + dev_k_Perl(i,1);
    k_HDPE(i+11,1) = k_plastic + dev_k_HDPE(i,1);
    h_Air(i+11,1)  = h + dev_h(i,1);
    L_Perl(i+11,1) = L_perlite + dev_L(i,1);
    L_HDPE(i+11,1) = L_plastic + dev_L(i,1);
end

k_Perl(11,1) = k_perlite;
k_HDPE(11,1) = k_plastic;
h_Air(11,1) = h;
L_Perl(11,1) = L_perlite;
L_HDPE(11,1) = L_plastic;

%% Analysis
% nom_index = [18;30]; % The nominal value index for first and second problem for the case study
TIME = [1800;3000];

[nominal_solution_fire,nominal_fire_Fo] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_plastic,k_plastic,h,k_perlite,L_perlite);
[nominal_solution_oven,nominal_oven_Fo] = HDE_solution(TIME(2,1),T_inf(2,1),T_s_i(2,1),L_plastic,k_plastic,h,k_perlite,L_perlite);

for i = 1:1:length(k_Perl)
    [solution_L_HDPE,L_HDPE_Fo,L_HDPE_Bi,L_HDPE_lam,L_HDPE_C1] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_HDPE(i,1),k_plastic,h,k_perlite,L_perlite);
    [solution_k_HDPE,k_HDPE_Fo,k_HDPE_Bi,k_HDPE_lam,k_HDPE_C1] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_plastic,k_HDPE(i,1),h,k_perlite,L_perlite);
    [solution_h_Air,h_Fo,h_Bi,h_lam,h_C1] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_plastic,k_plastic,h_Air(i,1),k_perlite,L_perlite);
    [solution_k_Perl,k_Perl_Fo,k_Perl_Bi,k_Perl_lam,k_Perl_C1] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_plastic,k_plastic,h,k_Perl(i,1),L_perlite);
    [solution_L_Perl,L_Perl_Fo,L_Perl_Bi,L_Perl_lam,L_Perl_C1] = HDE_solution(TIME(1,1),T_inf(1,1),T_s_i(1,1),L_plastic,k_plastic,h,k_perlite,L_Perl(i,1));
    solutions(i,1:5) = [solution_L_HDPE solution_k_HDPE solution_h_Air solution_k_Perl solution_L_Perl];
    sol_Fo(i,1:5) = [L_HDPE_Fo k_HDPE_Fo h_Fo k_Perl_Fo L_Perl_Fo];
    sol_Bi(i,1:5) = [L_HDPE_Bi k_HDPE_Bi h_Bi k_Perl_Bi L_Perl_Bi];
    sol_Lam(i,1:5)= [L_HDPE_lam k_HDPE_lam h_lam k_Perl_lam L_Perl_lam];
    sol_C1(i,1:5) = [L_HDPE_C1 k_HDPE_C1 h_C1 k_Perl_C1 L_Perl_C1];
   
    sensitivity_L_HDPE = nominal_solution_fire - solution_L_HDPE;
    sensitivity_k_HDPE = nominal_solution_fire - solution_k_HDPE;
    sensitivity_h_Air = nominal_solution_fire - solution_h_Air;
    sensitivity_k_Perl = nominal_solution_fire - solution_k_Perl;
    sensitivity_L_Perl = nominal_solution_fire - solution_L_Perl;
    sensitivities(i,1:5) = [sensitivity_L_HDPE sensitivity_k_HDPE sensitivity_h_Air sensitivity_k_Perl sensitivity_L_Perl];
end

percentage = [-100;-90;-80;-70;-60;-50;-40;-30;-20;-10;0;10;20;30;40;50;60;70;80;90;100];

%% Plotting
hold on
scatter(percentage,sensitivities(:,1),'d');
scatter(percentage,sensitivities(:,2),'s');
scatter(percentage,sensitivities(:,3),'x');
scatter(percentage,sensitivities(:,4),'o');
scatter(percentage,sensitivities(:,5),'+');
grid on
title('One-at-a-time Sensitivity Analysis');
ylabel('Temperature Difference from Nominal, \DeltaT [K]');
xlabel('Percent difference from Nominal Parameter Value');
legend('L_{HDPE}','k_{HDPE}','h_{Air}','k_{perlite}','L_{perlite}');

