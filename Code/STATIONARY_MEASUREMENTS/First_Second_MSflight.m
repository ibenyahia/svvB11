% Citation 550 - Linear simulation hai
% xcg = 0.25*c
load('MASS_CALCULATIONS\Total_mass_writtenvalues.mat')

% First Stationary flight
hp_1    = 0.3048*[9000,8990,8990,9000,9000,9010]; % pressure altitude [m]
V_CAS_1    = 0.514444*([250,224,190,164,130,119]-2);  % calibrated airspeed [m/s]
aoa_1 = [1.6,2.3,3.7,5.4,8.9,10.8]; % angle of attack [deg]
% th_1    = aoa_1; % pitch angle in the stationary flight condition [deg]
T_m_1 = 273.15 + [-0.8,-2.2,-4.7,-6.5,-7.9,-8.3]; % measured temperature [K]
FFl_1 = [724,629,495,452,430,427]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_1 = [766,674,534,486,475,465]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_1,p_1,V_EAS_1,V_TAS_1,delta_T_1,M_T_1] = getImpValues(hp_1,V_CAS_1,T_m_1); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]

% Second Stationary flight
hp_2    = 0.3048*[10960,11200,11380,11560,10850,10280,9800,10190,10220]; % pressure altitude [m]
V_CAS_2     = 0.514444*([157,145,136,127,164,173,183,157,157]-2);  % calibrated airspeed [m/s]
aoa_2 = [5.8,6.9,8,9.4,5.3,4.6,4,5.8,5.8]; % angle of attack [deg]
% th_2    = aoa_2; % pitch angle in the stationary flight condition [deg]
T_m_2 = 273.15 + [-10.4,-11.5,-12.5,-13.2,-9.8,-8.2,-7,-8.5,-8.8]; % measured temperature [K]
FFl_2 = [409,405,402,400,414,422,427,419,458]/3600*0.45359237; % fuel flow of the left engine [kg/s]
FFr_2 = [447,442,440,437,452,461,469,420,457]/3600*0.45359237; % fuel flow of the right engine [kg/s]
[rho_2,p_2,V_EAS_2,V_TAS_2,delta_T_2,M_T_2] = getImpValues(hp_2,V_CAS_2,T_m_2); % density [kg/m^3], pressure [Pa], equiv. airspeed [m/s], true airspeed [m/s], temp. diff. [K] and true mach [-]
de = [-0.4,-0.9,-1.4,-2,-0.1,0.2,0.5,-0.4,-0.8]; % Elevator deflection [deg]
Fe = [0,-24,-35,-47,17,39,60,1,-23]; % Elevator stick force [F]
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
% for i=1:7
%     fprintf("%f ",[hp_2(i),M_T_2(i),delta_T_2(i),FFl_2(i),FFr_2(i)]);
%     fprintf("\n")
% end
% Forces
Thr_L_1 = [3430.52,2929.69,2186.76,2021.7,2046.62,2088.28]; % Thrust of left engine [N]
Thr_R_1 = [3725.89,3252.46,2465.13,2273.13,2396.61,2388.78]; % Thrust of right engine [N]
Thr_tot_1 = Thr_L_1 + Thr_R_1; % Total thrust of the aircraft [N]

% Aerodynamic coefficients CL and CD
CD = Thr_tot_1./(0.5*rho_1.*V_TAS_1.^2*S); % Drag coefficient (T = D) [-]
CL = W(1:6)./(0.5*rho_1.*V_TAS_1.^2*S); % Lift coefficient (W = L) [-]
CL_sq = CL.^2; % Lift coefficient squared [-]

% -- Curve fit obtained values -- %

% CL^2 = CD*pi*A*e - CD0*pi*A*e: R^2 = 0.9961, RMSE = 0.02442
pAe = 19.22; % (95% confidence: 18.46 - 21.89) pi * A * e, obtained from curve fit: X_data = CD, Y_data = CL
e = pAe/(pi*A); % Oswald efficiency factor [-]
CD0 = 0.4109/(pAe); % (95% confidence: 0.3399 - 0.482) Zero-lift drag coefficient [-]

% CL = CLa*aoa + CL0 : R^2 = 0.9993, RMSE = 0.008558
CLa = 0.08011; % (95% confidence: 0.07723 - 0.08298) Lift coefficient gradient wrt angle of attack [deg^-1]
CL0 = 0.08326; % (95% confidence: 0.0628 - 0.09965) Lift coefficient at angle of attack = 0 [-]
aoa_0 = -1.01; % Angle of attack at zero lift [deg]
th0 = aoa_0*pi/180; % pitch angle
%%
%%% - - - Second Stationary Measurements - - - %%%

% Same calculations as in the second stationary measaurements

% Forces
Thr_L_2 = [1828.48,1872.66,1909.3,1951.76,1824.3,1810.54,1776.48]; % Thrust of left engine [N]
Thr_R_2 = [2111.71,2153.65,2201.34,2239.89,2104.44,2094.6,2078.57]; % Thrust of right engine [N]
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


% Longitudinal stability

Cmde = -0.5/(-0.8+0.4) *(W(14)+W(15))/(0.5*rho_2(7)*V_TAS_2(7)^2*S)*(x_cg(15)-x_cg(14))/c*180/pi; % Elevator effectiveness [rad^-1]
% Cmde = -1/(-0.5-0) * (CLa*5.3 + CL0)*(134-288)*0.0254/c; % Elevator effectiveness [deg^-1]
Cma = -Cmde*dde_da;% Moment coefficient slope wrt alpha [rad^-1]
%%

% for i=1:7
%     fprintf("%f ",[hp_2(i),M_T_2(i),delta_T_2(i),0.048,0.048]);
%     fprintf("\n")
% end
Ws = 60500; % Standard aircraft weight
V_EAS_r = V_EAS_2(1:7).*sqrt(Ws./W(7:13)); % Reduced Equivalent Airspeed
Thr_tot_r = 2*[1619.94,1690.32,1745.66,1801.23,1581.28,1514.02,1449.57]; % Reduced Thrust
eng_diam = 0.686; % Engine diameter from https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_JT15D
Thr_coeff = Thr_tot_2./(0.5*rho0*V_EAS_2(1:7).^2*pi*(eng_diam/2)^2); % Thrust coefficient
Thr_coeff_r = Thr_tot_r./(0.5*rho0*V_EAS_r.^2*pi*(eng_diam/2)^2); % Reduced Thrust coefficient
CmT = -0.0064; % Change in Cm wrt Tc (thrust coefficient)
de_red = de(1:7) - 1/Cmde*(Thr_coeff_r-Thr_coeff); % Reduced Elevator Deflection
Fe_red = Fe(1:7)./W(7:13) * Ws; % Reduced Elevator Control Force

%%

tmp_de = de(1:7);
tmp_de_red = de_red;
tmp_Fe = Fe(1:7);
tmp_Fe_red = Fe_red;
tmp_V_red = V_EAS_r;
[out,index] = sort(de_red);
for i = 1:length(out)
    tmp_de(i) = de(index(i));
    tmp_de_red(i) = de_red(index(i));
    tmp_Fe(i) = Fe(index(i));
    tmp_Fe_red(i) = Fe_red(index(i));
    tmp_V_red(i) = V_EAS_r(index(i));
end

figure(1)
hold on
sgtitle('Reduced Elevator Trim Curve (Reference Data)')
plot(tmp_V_red,tmp_de,['r','-o']);
plot(tmp_V_red,tmp_de_red,['b','-o']);
set(gca, 'YDir','reverse')
xlabel('$$\tilde{V}_{EAS}$$ [m/s]','Interpreter','Latex')
ylabel('\delta_e [rad]')
legend('\delta_{eq}','\delta_{eq}^{*}')
hold off

figure(2)
hold on
sgtitle('Reduced Elevator Control Force Curve (Reference Data)')
plot(tmp_V_red,tmp_Fe,['r','-o']);
plot(tmp_V_red,tmp_Fe_red,['b','-o']);
set(gca, 'YDir','reverse')
xlabel('$$\tilde{V}_{EAS}$$ [m/s]','Interpreter','Latex')
ylabel('F_{e} [N]')
legend('F_{e}','F_{e}^{*}')
hold off

%%

load census;
func_1 = fit(reshape(aoa_1,[length(aoa_1),1]),reshape(CL,[length(CL),1]),'poly1');
func_2 = fit(reshape(CD,[length(CD),1]),reshape(CL,[length(CL),1]),'poly2');
func_3 = fit(reshape(CD,[length(CD),1]),reshape(CL_sq,[length(CL_sq),1]),'poly1');
func_4 = fit(reshape(aoa_1,[length(aoa_1),1]),reshape(CD,[length(CD),1]),'poly2');


figure(3)

sgtitle('First Stationary Measurements (Flight Data)')
subplot(2,2,1)
plot(func_1,aoa_1,CL,'o');
title('C_{L} - \alpha')
xlabel('\alpha [deg]')
ylabel('C_{L} [-]')
legend('Polynomial of order 1','Data points')

subplot(2,2,2)
plot(func_2,CD,CL,'o');
title('Lift-Drag Polar C_{L} - C_{D}')
xlabel('C_{D} [-]')
ylabel('C_{L} [-]')
legend('Polynomial of order 2','Data points')

subplot(2,2,3)
plot(func_3,CD,CL_sq,'o')
title('C_{L}^{2} - C_{D}')
xlabel('C_{D} [-]')
ylabel('C_{L}^{2} [-]')
legend('Polynomial of order 1','Data points')

subplot(2,2,4)
plot(func_4,aoa_1,CD,'o');
title('C_{D} - \alpha')
xlabel('\alpha [deg]')
ylabel('C_{D} [-]')
legend('Polynomial of order 2','Data points')