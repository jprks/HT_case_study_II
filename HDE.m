function [T] = HDE(lam,C_1,T_i,T_inf,Fo)

T = T_inf + (T_i - T_inf)*C_1*exp(-lam^2*Fo);

end