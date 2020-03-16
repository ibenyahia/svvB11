% Citation 550 - Linear simulation hai

% xcg = 0.25*c

% First Stationary flight
hp_1    = 0.3048*[5010,5020,5020,5030,5020,5110]; % pressure altitude [m]
V_CAS_1    = 0.514444*([249,221,192,163,130,118]-2);  % calibrated airspeed [m/s]
aoa_1 = [1.7,2.4,3.6,5.4,8.7,10.6]; % angle of attack [deg]
th_1    = aoa_1; % pitch angle in the stationary flight condition [deg]
T_m_1 = 273.15 + [12.5,10.5,8.8,7.2,6,5.2]; % measured temperature [K]
FFl_1 = [798,673,561,463,443,474]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_1 = [813,682,579,484,467,499]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_1,p_1,V_EAS_1,V_TAS_1,delta_T_1,M_T_1] = getImpValues(hp_1,V_CAS_1,T_m_1); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]

% Second Stationary flight
hp_2    = 0.3048*[6060,6350,6550,6880,6160,5810,5310]; % pressure altitude [m]
V_CAS_2     = 0.514444*([161,150,140,130,173,179,192]-2);  % calibrated airspeed [m/s]
aoa_2 = [5.3,6.3,7.3,8.5,4.5,4.1,3.4]; % angle of attack [deg]
th_2    = aoa_2; % pitch angle in the stationary flight condition [deg]
T_m_2 = 273.15 + [5.5,4.5,3.5,2.5,5,6.2,8.2]; % measured temperature [K]
FFl_2 = [462,458,454,449,465,472,482]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_2 = [486,482,477,473,489,496,508]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_2,p_2,V_EAS_2,V_TAS_2,delta_T_2,M_T_2] = getImpValues(hp_2,V_CAS_2,T_m_2); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]
de = [0,-0.4,-0.9,-1.5,0.4,0.6,1]; % Elevator deflection [deg]
Fe = [0,-23,-29,-46,26,40,83]; % Elevator stick force [F]
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
Thr_L_1 = [3609.17,2961.13,2379.6,1852.96,1888.86,2207.36]; % Thrust of left engine [N]
Thr_R_1 = [3711.33,3022.38,2505.27,2004.88,2072.03,2403.69]; % Thrust of right engine [N]
Thr_tot_1 = Thr_L_1 + Thr_R_1; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
CD_1 = Thr_tot_1./(0.5*rho0*V_EAS_1.^2*S); % Drag coefficient (T = D) [-]
CL_1 = W(1:6)./(0.5*rho0*V_EAS_1.^2*S); % Lift coefficient (W = L) [-]
CL_sq_1 = CL_1.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9807, RMSE = 0.05639
pAe_1 = 19.67; % (95% confidence: 15.84 - 23.51) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
e_1 = pAe_1/(pi*A); % Oswald efficiency factor [-]
CD0_1 = 0.4107/(pAe_1); % (95% confidence: 0.2486 - 0.5728) Zero-lift drag coefficient [-]

% CL = CLa*aoa + CL0 : R^2 = 0.999, RMSE = 0.01039
CLa_1 = 0.084; % (95% confidence: 0.08039 - 0.08761) Lift coefficient gradient wrt angle of attack [deg^-1]
CL0_1 = 0.0738; % (95% confidence: 0.05103 - 0.09657) Lift coefficient at angle of attack = 0 [-]
aoa_0_1 =  -0.8726; % (95% confidence: -1.176 - -0.5688) Angle of attack at zero lift [deg]
%%
%%% - - - Second Stationary Measurements - - - %%%

% Same calculations as in the second stationary measaurements

% Forces
Thr_L_2 = [1907.6,1954.58,1993.6,2033.12,1871.72,1870.59,1844.14]; % Thrust of left engine [N]
Thr_R_2 = [2082.39,2132.35,2166.83,2214.97,2042.99,2039.63,2024.19]; % Thrust of right engine [N]
Thr_tot_2 = Thr_L_2 + Thr_R_2; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
CD_2 = Thr_tot_2./(0.5*rho0*V_EAS_2.^2*S); % Drag coefficient (T = D) [-]
CL_2 = W(7:13)./(0.5*rho0*V_EAS_2.^2*S); % Lift coefficient (W = L) [-]
CL_sq_2 = CL_2.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9946, RMSE = 0.01436
pAe_2 = 15.73; % (95% confidence: 14.4 - 17.06) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
e_2 = pAe_2/(pi*A); % Oswald efficiency factor [-]
CD0_2 = 0.2413/(pAe_2); % (95% confidence: 0.1921 - 0.2906) Zero-lift drag coefficient [-]

% CL = CLa*(aoa - aoa_0 = 0) + CL0: R^2 = 0.9994, RMSE = 0.004122
CLa_2 = 0.08431; % (95% confidence: 0.08195 - 0.08667) Lift coefficient gradient wrt angle of attack [deg^-1]
CL0_2 = 0.06981; % (95% confidence: 0.05594 - 0.08368) Lift coefficient at angle of attack = 0 [-]
aoa_0_2 =  -0.8242; % (95% confidence: -1.011 - -0.6375) Angle of attack at zero lift [deg]

% R^2 = 0.9986, RMSE = 0.03652
dde_da = -0.4797; % (95% confidence: -0.5006 - -0.4588) Slope elevator deflection wrt to angle of attack

% Longitudinal stability

Cmde = -1/(-0.5-0) * (CLa_2*5.3 + CL0_2)*(134-288)*0.0254/c; % Elevator effectiveness [deg^-1]
Cma = -Cmde*dde_da; % Moment coefficient slope wrt alpha [deg^-1]
%%

V_red_1 = V_EAS_1*sqrt();

 
%%
% Constant values concerning aircraft inertia

muc_1  = m./(rho_1.*S*c); % (FSM)
muc_2  = m./(rho_2.*S*c); % (SSM)
mub_1  = m./(rho_1.*S*b); % (FSM)
mub_2  = m./(rho_2.*S*b); % (SSM)
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa_1 = CLa_1;   		        % Wing normal force slope [ ] (FSM)
CNwa_2 = CLa_2;                 % Wing normal force slope [ ] (SSM)
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

% CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
% CD = CD0 + (CLa*aoa_1)^2/(pi*A*e);  % Drag coefficient [ ]

% Stability derivatives

CX0_1  = W(1:6).*sin(th_1)/(0.5*rho_1.*V_EAS_1.^2*S); % (FSM)
CX0_2  = W(7:13).*sin(th_2)/(0.5*rho_2.*V_EAS_2.^2*S); % (SSM)
CXu    = -0.02792;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0_1  = -W(1:6).*cos(th_1)/(0.5*rho_1.*V_EAS_1.^2*S);  % (FSM)
CZ0_2  = -W(7:13).*cos(th_2)/(0.5*rho_2.*V_EAS_2.^2*S); % (SSM)
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

% As = array([]);
