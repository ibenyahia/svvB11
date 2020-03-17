% Citation 550 - Linear simulation hai

% xcg = 0.25*c

% First Stationary flight
hp_1    = 0.3048*[5010,5020,5020,5030,5020,5110]; % pressure altitude [m]
V_CAS_1    = 0.514444*([249,221,192,163,130,118]-2);  % calibrated airspeed [m/s]
aoa_1 = [1.7,2.4,3.6,5.4,8.7,10.6]; % angle of attack [deg]
% th_1    = aoa_1; % pitch angle in the stationary flight condition [deg]
T_m_1 = 273.15 + [12.5,10.5,8.8,7.2,6,5.2]; % measured temperature [K]
FFl_1 = [798,673,561,463,443,474]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_1 = [813,682,579,484,467,499]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_1,p_1,V_EAS_1,V_TAS_1,delta_T_1,M_T_1] = getImpValues(hp_1,V_CAS_1,T_m_1); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]

% Second Stationary flight
hp_2    = 0.3048*[6060,6350,6550,6880,6160,5810,5310,5730,5790]; % pressure altitude [m]
V_CAS_2     = 0.514444*([161,150,140,130,173,179,192,161,161]-2);  % calibrated airspeed [m/s]
aoa_2 = [5.3,6.3,7.3,8.5,4.5,4.1,3.4,5.3,5.3]; % angle of attack [deg]
% th_2    = aoa_2; % pitch angle in the stationary flight condition [deg]
T_m_2 = 273.15 + [5.5,4.5,3.5,2.5,5,6.2,8.2,5,5]; % measured temperature [K]
FFl_2 = [462,458,454,449,465,472,482,471,493]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_2 = [486,482,477,473,489,496,508,468,490]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_2,p_2,V_EAS_2,V_TAS_2,delta_T_2,M_T_2] = getImpValues(hp_2,V_CAS_2,T_m_2); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]
de = [0,-0.4,-0.9,-1.5,0.4,0.6,1,0,-0.5]; % Elevator deflection [deg]
Fe = [0,-23,-29,-46,26,40,83,0,-30]; % Elevator stick force [F]
detr = 2.8; % Elevator trim [deg]

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
Thr_L_1 = [3586.64,2947.83,2372.32,1849.75,1887.69,2206.97]; % Thrust of left engine [N]
Thr_R_1 = [3688.35,3008.91,2497.77,2001.51,2070.79,2403.28]; % Thrust of right engine [N]
Thr_tot_1 = Thr_L_1 + Thr_R_1; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
CD = Thr_tot_1./(0.5*rho_1.*V_TAS_1.^2*S); % Drag coefficient (T = D) [-]
CL = W(1:6)./(0.5*rho_1.*V_TAS_1.^2*S); % Lift coefficient (W = L) [-]
CL_sq = CL.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9826, RMSE = 0.05349
pAe = 19.52; % (95% confidence: 15.91 - 23.12) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
e = pAe/(pi*A); % Oswald efficiency factor [-]
CD0 = 0.4015/(pAe); % (95% confidence: 0.25 - 0.553) Zero-lift drag coefficient [-]

% CL = CLa*aoa + CL0 : R^2 = 0.999, RMSE = 0.01052
CLa = 0.08446; % (95% confidence: 0.0808 - 0.08811) Lift coefficient gradient wrt angle of attack [deg^-1]
CL0 = 0.06751; % (95% confidence: 0.04445 - 0.09057) Lift coefficient at angle of attack = 0 [-]
aoa_0 = -0.7933; % Angle of attack at zero lift [deg]
th0 = aoa_0*pi/180; % pitch angle


%%
%%% - - - Second Stationary Measurements - - - %%%

% Same calculations as in the second stationary measaurements

% Forces
Thr_L_2 = [1903.82,1951.57,1991.56,2031.63,1868.25,1866.2,1837.87]; % Thrust of left engine [N]
Thr_R_2 = [2078.41,2129.19,2164.49,2213.41,2039.33,2035.04,2017.57]; % Thrust of right engine [N]
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

% R^2 = 0.9986, RMSE = 0.03652
dde_da = -0.4797; % (95% confidence: -0.5006 - -0.4588) Slope elevator deflection wrt to angle of attack

% Longitudinal stability

Cmde = -0.5/(-0.5) *(W(14)+W(15))/(0.5*rho0*V_EAS_2(8)^2*S)*(x_cg(15)-x_cg(14))/c*180/pi ; % Elevator effectiveness [rad^-1]
% Cmde = -1/(-0.5-0) * (CLa*5.3 + CL0)*(x_cg(15)-x_cg(14))*0.0254/c*180/pi; % Elevator effectiveness [deg^-1]
Cma = -Cmde*dde_da; % Moment coefficient slope wrt alpha [deg^-1]
%%

Ws = 60500; % Standard aircraft weight
V_EAS_r = V_EAS_2(1:7).*sqrt(Ws./W(7:13)); % Reduced Equivalent Airspeed
Thr_tot_r = 2*[1331.19,1394.62,1452.12,1516.85,1279.09,1242.83,1165.44]; % Reduced Thrust
eng_diam = 0.686; % Engine diameter from https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_JT15D
Thr_coeff = Thr_tot_2./(0.5*rho0*V_EAS_2(1:7).^2*pi*(eng_diam/2)^2); % Thrust coefficient
Thr_coeff_r = Thr_tot_r./(0.5*rho0*V_EAS_r.^2*pi*(eng_diam/2)^2); % Reduced Thrust coefficient
CmT = -0.0064; % Change in Cm wrt Tc (thrust coefficient)
de_red = de(1:7) - 1/Cmde*(Thr_coeff_r-Thr_coeff); % Reduced Elevator Deflection
Fe_red = Fe(1:7)./W(7:13) * Ws; % Reduced Elevator Control Force