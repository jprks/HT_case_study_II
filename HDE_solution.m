function [nominal_temp,nominal_Fo,nominal_Bi,nominal_Lam,nominal_C1] = HDE_solution(t,T_inf,T_i,L_PLASTIC,K_PLASTIC,H_AIR,K_PERLITE,L_PERLITE)
global cp rho
load coefficients
R_eff_prime = 1/H_AIR + L_PLASTIC/K_PLASTIC; % [K/W] - Effective thermal resistance of convection and conduction on plastic wall
U = 1/R_eff_prime; % [W/K] - Transforming R_eff_prime into a form that can be used for the HDE solution
alp = K_PERLITE/(cp*rho); % [m^2/s] - Thermal diffusivity
L_c = L_PERLITE+L_PLASTIC; % [m] - Characteristic length
Bi = U*L_c/K_PERLITE; % [-] - Biot number

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

Fo = alp*t/L_c^2;

Temp = HDE(lam,C_1,T_i,T_inf,Fo);

nominal_temp = Temp;
nominal_Fo = Fo;
nominal_Bi = Bi;
nominal_Lam = lam;
nominal_C1 = C_1;

end