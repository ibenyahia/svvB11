% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

% Get variables 
t = flightdata.time.data;                        % time [sec]
hp = flightdata.Dadc1_alt.data*0.3048;           % altitude [m]
alpha = flightdata.vane_AOA.data*pi/180 ;        % angle of attack [rad]
p = flightdata.Ahrs1_bRollRate.data*pi/180;      % roll rate [rad/s]
q = flightdata.Ahrs1_bPitchRate.data*pi/180;     % pitch rate [rad/s]
r = flightdata.Ahrs1_bYawRate.data*pi/180;       % yaw rate [rad/s]
phi = flightdata.Ahrs1_Roll.data*pi/180;         % roll [rad]
theta = flightdata.Ahrs1_Pitch.data*pi/180 ;     % pitch [rad]
delta_a = flightdata.delta_a.data*pi/180;        % aileron deflection [rad]
delta_e = flightdata.delta_e.data*pi/180 ;       % elevator deflection [rad]
delta_r = flightdata.delta_r.data*pi/180;        % rudder deflection [rad]
vtas = flightdata.Dadc1_tas.data*0.514444;       % true airspeed [m/s]  

% Index for reference data
i_SPref = 36241;
i_PHref = 32201; % 32281
i_DRref = 37081;
i_DRDampref = 37571;
i_APref = 35411;
i_SPIRref = 39111;

% Index for flight data

i_SPfd = 30401;
i_PHfd = 32511; %32511
i_DRfd = 34411;
i_DRDampfd = 35211;
i_APfd = 31611;
i_SPIRfd = 36781; % 36511

index = i_PHfd; % Index of interest for initial conditions

%Order in list: hp0 [m], V0[m/s], alpha0 [rad], th0 [rad], mass [kg],td [s]

% Reference data
shortperiod_ref = {hp(i_SPref),vtas(i_SPref),alpha(i_SPref),theta(i_SPref),Mtotal(i_SPref),8};
phugoid_ref = {hp(i_PHref),vtas(i_PHref),alpha(i_PHref),theta(i_PHref),Mtotal(i_PHref),228};
dr_ref = {hp(i_DRref),vtas(i_DRref),alpha(i_DRref),theta(i_DRref),Mtotal(i_DRref),19};
drdamp_ref = {hp(i_DRDampref),vtas(i_DRDampref),alpha(i_DRDampref),theta(i_DRDampref),Mtotal(i_DRDampref),11};
apr_ref = {hp(i_APref),vtas(i_APref),alpha(i_APref),theta(i_APref),Mtotal(i_APref),12};
spiral_ref = {hp(i_SPIRref),vtas(i_SPIRref),alpha(i_SPIRref),theta(i_SPIRref),Mtotal(i_SPIRref),200};

%Flight data
shortperiod_fd = {hp(i_SPfd),vtas(i_SPfd),alpha(i_SPfd),theta(i_SPfd),Mtotal(i_SPfd),6};
phugoid_fd = {hp(i_PHfd),vtas(i_PHfd),alpha(i_PHfd),theta(i_PHfd),Mtotal(i_PHfd), 150};
dr_fd = {hp(i_DRfd),vtas(i_DRfd),alpha(i_DRfd),theta(i_DRfd),Mtotal(i_DRfd), 25};
drdamp_fd = {hp(i_DRDampfd),vtas(i_DRDampfd),alpha(i_DRDampfd),theta(i_DRDampfd),Mtotal(i_DRDampfd),11};
apr_fd = {hp(i_APfd),vtas(i_APfd),alpha(i_APfd),theta(i_APfd),Mtotal(i_APfd), 14};
spiral_fd = {hp(i_SPIRfd),vtas(i_SPIRfd),alpha(i_SPIRfd),theta(i_SPIRfd),Mtotal(i_SPIRfd) , 145};

selection = phugoid_fd;   %Replace name with flight condition of interest

% Initial conditions
hp0    = selection{1}(1);             % Initial height [m]
V0     = selection{2}(1);             % Initial airspeed [m/s]
alpha0 = 0;                           % Alpha stability-axis [rad]
th0    = 0;                           % Theta stability-axis [rad]
m      = selection{5}(1);         	  % Initial aircraft mass [kg]
td     = selection{6};                % Eigenmotion duration [s]
% aerodynamic properties

% order: Oswald factor [-], Zero lift drag coefficient [-], CLa [rad^-1]
aerocoeff_ref = {0.7314,0.0208,4.8111};  % Reference data
aerocoeff_flight = {0.725,0.0214,4.59};  % Flight data
aerocoeff_select = aerocoeff_ref;     % Select data of interest

e      = aerocoeff_select{1};              % Oswald factor [-]
CD0    = aerocoeff_select{2};             % Zero lift drag coefficient [-]
CLa    = aerocoeff_select{3};            % Slope of CL-alpha curve [rad^-1]

% Longitudinal stability
coef_fd = {-0.5347, -1.1494};       % Flight data
coef_ref = {-0.5718, -1.1935};      % Reference data
coef_select = coef_ref;              % Select data of interest

Cma    = coef_select{1};            % Longitudinal stability [ ]
Cmde   = coef_select{2};            % Elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % Wing area [m^2]
Sh     = 0.2*S;           % Stabiliser area [m^2]
Sh_S   = Sh/S;	          % [-]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [-]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
V0h_V0   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant V0alues concerning atmosphere and graV0ity

rho0   = 1.2250;          % air density at sea leV0el [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea leV0el in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (graV0ity constant)


rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

% Constant Values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.3925;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CL)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W.*sin(th0)/(0.5*rho.*V0.^2.*S);
CXu    = -0.095;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W.*cos(th0)/(0.5*rho.*V0.^2.*S);
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


%Numerical model stuff starts here
%V0ARIABLES THAT ARE MISSING

CXdt = 0.000;   %trim tab is really small for this aircraft so it's contribution
CZdt = 0.000;   %compared to the actual elevator is an order of 
Cmdt = 0.000;   %magnitude smaller

%%

%%%% Entries of symmetric eigenmotions

as_11 = -2*muc*c/(V0^2);
as_22 = (CZadot-2*muc)*c/V0;
as_33 = -c/V0;
as_42 = Cmadot*c/V0;
as_44 = -2*muc*KY2*(c/V0)^2;

bs_11 = CXu/V0;
bs_12 = CXa;
bs_13 = CZ0;
bs_14 = CXq*c/V0;
bs_21 = CZu/V0;
bs_22 = CZa;
bs_23 = -CX0;
bs_24 = (CZq+2*muc)*c/V0;
bs_34 = c/V0;
bs_41 = Cmu/V0;
bs_42 = Cma;
bs_44 = Cmq*c/V0;

cs_11 = CXde;
cs_21 = CZde;
cs_41 = Cmde;

C1_s = [as_11 0 0 0; 0 as_22 0 0; 0 0 as_33 0; 0 as_42 0 as_44];
C2_s = [bs_11 bs_12 bs_13 bs_14;bs_21 bs_22 bs_23 bs_24; 0 0 0 bs_34; bs_41 bs_42 0 bs_44];
C3_s = [cs_11;cs_21;0;cs_41];

% Symmetric state space system
A_s = -1*C1_s\C2_s;
B_s = -1*C1_s\C3_s;
C_s = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% D_s = zeros(4,1);
D_s = -1.5*ones(4,1);
% D_s = ones(4,1);
sym_sys = ss(A_s,B_s,C_s,D_s); % State-space system
eig_sym = eig(A_s); % Eigenvalues


%%%% Entries of asymmetric eigenmotions
aa_11 = (CYbdot-2*mub)*b/V0;
aa_22 = -b/(2*V0);
aa_33 = -2*mub*KX2*(b/V0)^2;
aa_34 = 2*mub*KXZ*(b/V0)^2;
aa_41 = Cnbdot*b/V0;
aa_43 = 2*mub*KXZ*(b/V0)^2;
aa_44 = -2*mub*KZ2*(b/V0)^2;

ba_11 = CYb;
ba_12 = CL;
ba_13 = CYp*b/(2*V0);
ba_14 = (CYr-4*mub)*b/(2*V0);
ba_23 = b/(2*V0);
ba_31 = Clb;
ba_33 = Clp*b/(2*V0);
ba_34 = Clr*b/(2*V0);
ba_41 = Cnb;
ba_43 = Cnp*b/(V0*2);
ba_44 = Cnr*b/(V0*2);

ca_11 = CYda;
ca_12 = CYdr;
ca_31 = Clda;
ca_32 = Cldr;
ca_41 = Cnda;
ca_42 = Cndr;

C1_a = [aa_11 0 0 0; 0 aa_22 0 0; 0 0 aa_33 aa_34; aa_41 0 aa_43 aa_44];
C2_a = [ba_11 ba_12 ba_13 ba_14; 0 0 ba_23 0; ba_31 0 ba_33 ba_34; ba_41 0 ba_43 ba_44];
C3_a = [ca_11 ca_12; 0 0; ca_31 ca_32; ca_41 ca_42];

% Asymmetric state space system
A_a = -inv(C1_a)*C2_a;
B_a = -inv(C1_a)*C3_a;
C_a = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
D_a = zeros(4,2);

asym_sys = ss(A_a,-B_a,C_a,D_a); % State-space
eig_asym = eig(A_a); % Eigenvalues

%%
time = linspace(0,td,td*10+1);
endind = index+length(time)-1;

sym_resp = lsim(sym_sys,delta_e(index:endind),time,[0,0,0,q(index)]);
vtas_resp = sym_resp(:,1);
alpha_resp = sym_resp(:,2);
theta_resp = sym_resp(:,3);
q_resp = sym_resp(:,4);


alphastab = alpha(index:endind)-alpha(index);
thetastab = theta(index:endind)-theta(index);
vtasstab = vtas(index:endind)-vtas(index);


asym_resp = lsim(asym_sys,[delta_a(index:endind),delta_r(index:endind)],time,[0,phi(index),p(index),r(index)]);
phi_resp = asym_resp(:,2);
p_resp = asym_resp(:,3);
r_resp = asym_resp(:,4);

% short_period_plot(alphastab,thetastab,delta_e(index:endind),q(index:endind),time,alpha_resp,q_resp,time,theta_resp);
% phugoid_plot(vtas(index:endind),thetastab,delta_e(index:endind),q(index:endind),time,vtas_resp+V0,theta_resp,q_resp,time);
% dr_plot(p(index:endind),phi(index:endind),r(index:endind),delta_a(index:endind),delta_r(index:endind),time,p_resp,phi_resp,r_resp,time,1);
% spiral_plot(p(index:endind),phi(index:endind),r(index:endind),delta_a(index:endind),delta_r(index:endind),time,p_resp,phi_resp,r_resp,time);
% ap_roll_plot(p(index:endind),phi(index:endind),r(index:endind),delta_a(index:endind),delta_r(index:endind),time,p_resp,phi_resp,r_resp,time);

%%

time_2 = linspace(0,150,1500);
cont_resp = lsim(sym_sys,-10*pi/180*ones(length(time_2),1),time_2,[0,0,0,0]);

figure(2)
sgtitle(' Control Input Response \delta_e = -1.5 \circ')

subplot(4,2,1);
plot(time_2,cont_resp(:,1)+V0)
title('True Airspeed V_{TAS} versus Time t')
ylabel('V_{TAS} [m/s]')
xlabel('t [s]')

subplot(4,2,3);
plot(time_2,cont_resp(:,2))
title('Angle of Attack \alpha versus Time t')
ylabel('\alpha [rad]')
xlabel('t [s]')

subplot(4,2,5);
plot(time_2,cont_resp(:,3))
title('Pitch Angle \theta versus Time t')
ylabel('\theta [rad]')
xlabel('t [s]')

subplot(4,2,7);
plot(time_2,cont_resp(:,4))
title('Pitch Rate q versus Time t')
ylabel('q [rad/s]')
xlabel('t [s]')

subplot(4,2,[2,4])
plot(time_2,cont_resp(:,2))
ylabel('\alpha [rad]')
xlabel('t [s]')
xlim([0,6])

subplot(4,2,[6,8])
plot(time_2,cont_resp(:,4))
ylabel('q [rad/s]')
xlabel('t [s]')
xlim([0,6])



