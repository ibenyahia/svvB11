% Citation 550 - Linear simulation hai

% xcg = 0.25*c

% First Stationary flight
hp_1    = 0.3048*[9000,8990,8990,9000,9000,9010]; % pressure altitude [m]
V_CAS_1    = 0.514444*([250,224,190,164,130,119]-2);  % calibrated airspeed [m/s]
aoa_1 = [1.6,2.3,3.7,5.4,8.9,10.8]; % angle of attack [deg]
% th_1    = aoa_1; % pitch angle in the stationary flight condition [deg]
T_m_1 = 273.15 + [-0.8,-2.2,-4.7,-6.5,-7.9,-8.3]; % measured temperature [K]
FFl_1 = [724,629,495,452,430,427]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_1 = [766,674,534,486,475,465]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_1,p_1,V_EAS_1,V_TAS_1,delta_T_1,M_T_1] = getImpValues(hp_1,V_CAS_1,T_m_1); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]

% for i=1:6
%     fprintf("%f ",[hp_1(i),M_T_1(i),delta_T_1(i),FFl_1(i),FFr_1(i)]);
%     fprintf("\n")
% end

% Second Stationary flight
hp_2    = 0.3048*[10960,11200,11380,11560,10850,10280,9800]; % pressure altitude [m]
V_CAS_2     = 0.514444*([157,145,136,127,164,173,183]-2);  % calibrated airspeed [m/s]
aoa_2 = [5.8,6.9,8,9.4,5.3,4.6,4]; % angle of attack [deg]
% th_2    = aoa_2; % pitch angle in the stationary flight condition [deg]
T_m_2 = 273.15 + [-10.4,-11.5,-12.5,-13.2,-9.8,-8.2,-7]; % measured temperature [K]
FFl_2 = [409,405,402,400,414,422,427]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_2 = [447,442,440,437,452,461,469]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_2,p_2,V_EAS_2,V_TAS_2,delta_T_2,M_T_2] = getImpValues(hp_2,V_CAS_2,T_m_2); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]
de = [-0.4,-0.9,-1.4,-2,-0.1,0.2,0.5]; % Elevator deflection [deg]
Fe = [0,-24,-35,-47,17,39,60]; % Elevator stick force [F]
detr = 2.5; % Elevator trim [deg]

% Aircraft mass
m      = Mtotal; % mass [kg]

%%

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
%%

W = reshape(m*g,1,length(m)); % Aircraft weight [N]

%%% - - - First Stationary Measurements - - - %%%

% Forces
Thr_L_1 = [3345.18,2861.91,2144.25,1987.92,2021.91,2066.12]; % Thrust of left engine [N]
Thr_R_1 = [3634.85,3182.32,2419.77,2236.78,2369.37,2367.77]; % Thrust of right engine [N]
Thr_tot_1 = Thr_L_1 + Thr_R_1; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
CD = Thr_tot_1./(0.5*rho_1.*V_TAS_1.^2*S); % Drag coefficient (T = D) [-]
CL = W(1:6)./(0.5*rho_1.*V_TAS_1.^2*S); % Lift coefficient (W = L) [-]
CL_sq = CL.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9963, RMSE = 0.02569
pAe = 20.17; % (95% confidence: 18.46 - 21.89) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
e = pAe/(pi*A); % Oswald efficiency factor [-]
CD0 = 0.4339/(pAe); % (95% confidence: 0.3599 - 0.5078) Zero-lift drag coefficient [-]

% CL = CLa*aoa + CL0 : R^2 = 0.9993, RMSE = 0.009267
CLa = 0.08323; % (95% confidence: 0.08011 - 0.08634) Lift coefficient gradient wrt angle of attack [deg^-1]
CL0 = 0.08326; % (95% confidence: 0.06331 - 0.1032) Lift coefficient at angle of attack = 0 [-]
aoa_0 = -0.9957; % Angle of attack at zero lift [deg]
th0 = aoa_0*pi/180; % pitch angle
%%
%%% - - - Second Stationary Measurements - - - %%%

% Same calculations as in the second stationary measaurements

% Forces
Thr_L_2 = [1797.6,1844.57,1882.79,1927.23,1790.95,1776.43,1740]; % Thrust of left engine [N]
Thr_R_2 = [2078.22,2122.99,2172.38,2215.22,2069.1,2058.12,2038.24]; % Thrust of right engine [N]
Thr_tot_2 = Thr_L_2 + Thr_R_2; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
% CD_2 = Thr_tot_2./(0.5*rho0*V_EAS_2.^2*S); % Drag coefficient (T = D) [-]
% CL_2 = W(7:13)./(0.5*rho0*V_EAS_2.^2*S); % Lift coefficient (W = L) [-]
% CL_sq_2 = CL_2.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9946, RMSE = 0.01436
% pAe_2 = 15.73; % (95% confidence: 14.4 - 17.06) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
% e_2 = pAe_2/(pi*A); % Oswald efficiency factor [-]
% CD0_2 = 0.2413/(pAe_2); % (95% confidence: 0.1921 - 0.2906) Zero-lift drag coefficient [-]

% CL = CLa*(aoa - aoa_0 = 0) + CL0: R^2 = 0.9994, RMSE = 0.004122
% CLa_2 = 0.08431; % (95% confidence: 0.08195 - 0.08667) Lift coefficient gradient wrt angle of attack [deg^-1]
% CL0_2 = 0.06981; % (95% confidence: 0.05594 - 0.08368) Lift coefficient at angle of attack = 0 [-]
% aoa_0_2 =  -0.8242; % (95% confidence: -1.011 - -0.6375) Angle of attack at zero lift [deg]
% th0_2 = a0a_0_2*np.pi/180; % pitch angle

% R^2 = 0.9989, RMSE = 0.03208
dde_da = -0.4652; % (95% confidence: -0.4827 - -0.4477) Slope elevator deflection wrt to angle of attack

% for i=1:7
%     fprintf("%f ",[hp_2(i),M_T_2(i),delta_T_2(i),FFl_2(i),FFr_2(i)]);
%     fprintf("\n")
% end
% Longitudinal stability

Cmde = -0.5/(-0.8+0.4) *(W(14)+W(15))/(0.5*rho_2(7)*V_TAS_2(7)^2*S)*(x_cg(15)-x_cg(14))/c*180/pi ; % Elevator effectiveness [deg^-1]
% Cmde = -1/(-0.5-0) * (CLa*5.3 + CL0)*(134-288)*0.0254/c; % Elevator effectiveness [deg^-1]
Cma = -Cmde*dde_da; % Moment coefficient slope wrt alpha [deg^-1]
%%

% for i=1:7
%     fprintf("%f ",[hp_2(i),M_T_2(i),delta_T_2(i),0.048,0.048]);
%     fprintf("\n")
% end

Ws = 60500; % Standard aircraft weight
V_EAS_r = V_EAS_2.*sqrt(Ws./W(7:13)); % Reduced Equivalent Airspeed
Thr_tot_r = 2*[1591.62,1664,1720.88,1777.91,1551.59,1483.58,1417.51]; % Reduced Thrust
eng_diam = 0.686; % Engine diameter from https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_JT15D
Thr_coeff = Thr_tot_2./(0.5*rho0*V_EAS_2.^2*pi*(eng_diam/2)^2); % Thrust coefficient
Thr_coeff_r = Thr_tot_r./(0.5*rho0*V_EAS_r.^2*pi*(eng_diam/2)^2); % Reduced Thrust coefficient
CmT = -0.0064; % Change in Cm wrt Tc (thrust coefficient)
de_red = de - 1/Cmde*(Thr_coeff_r-Thr_coeff); % Reduced Elevator Deflection
Fe_red = Fe./W(7:13) * Ws; % Reduced Elevator Control Force


figure(1)
plot(V_EAS_r,de_red,'b');
hold on
plot(V_EAS_r,de,'r');
set(gca, 'YDir','reverse')
hold off


