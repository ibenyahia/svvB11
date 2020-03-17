function [rho,p,V_EAS,V_TAS,delta_T,M_T] = getImpValues(h_p,V_CAS,T_meas)
    
% constants parameters:
    R_c = 287.05;
    g = 9.80665;
    T_0 = 288.15;
    gamma = 1.4;
    rho_0 = 1.225;
    p_0 = 101325;
    lam = -0.0065;
% calculations
    T_ISA = T_0 + lam*h_p;
    rho = rho_0*((1+(lam*h_p/T_0))).^(-((g/(lam*R_c))+1));
    p = p_0*(1+lam*h_p/T_0).^(-g/(lam*R_c));
    M_T = zeros(1,length(p));
    for i=1:length(p)
    M_T(:,i) = sqrt(2/(gamma-1)*((1+p_0/p(i)*((1+(gamma-1)/(2*gamma) * rho_0/p_0 * V_CAS(i)^2)^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1));
    end
    T_true = T_meas./(1+(gamma-1)/2 * M_T.^2);
    a_speed = sqrt(gamma*R_c*T_true);
    V_TAS = M_T.*a_speed;
    V_EAS = V_TAS.*sqrt(rho/rho_0);
    delta_T = T_meas - T_ISA;
end


