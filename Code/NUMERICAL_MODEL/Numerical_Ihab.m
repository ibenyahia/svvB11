% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

%get variables 
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
i_PHref = 32281;
i_DRref = 37081;
i_DRDampref = 37571;
i_APref = 35411;
i_SPIRref = 39111;

% Index for flight data
i_SPfd = 30401;
i_PHfd = 32511;
i_DRfd = 34411;
i_DRDampfd = 35211;
i_APfd = 31611;
i_SPIRfd = 36781; % 36511

index = i_SPIRfd;

%Order in list: hp0 [m], V0[m/s], alpha0 [rad], th0 [rad], mass [kg],td [s]
shortperiod_ref = {hp(i_SPref),vtas(i_SPref),alpha(i_SPref),theta(i_SPref),Mtotal(i_SPref),8};
phugoid_ref = {hp(i_PHref),vtas(i_PHref),alpha(i_PHref),theta(i_PHref),Mtotal(i_PHref),220};
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

%Order in list: hp0 [m], V0[m/s], alpha0 [rad], th0 [rad], mass [kg], INDEX
selection = shortperiod_fd;   %Replace name with flight condition of interest

hp0    = selection{1}(0);  	  
V0     = selection{2}(0);    
alpha0 = selection{3}(0);        	  
th0    = selection{4}(0);   

% Aircraft mass
m      = selection{5}(0);         	  % mass [kg]

% aerodynamic properties
% order: Oswald factor [-], Zero lift drag coefficient [-], CLa [rad^-1]
aerocoeff_ref = {0.7314,0.0208,4.8111};
aerocoeff_flight = {0.725,0.0214,4.59};
aerocoeff_select = aerocoeff_ref;

e      = aerocoeff_select{1};              % Oswald factor [ ]
CD0    = aerocoeff_select{2};             % Zero lift drag coefficient [ ]
CLa    = aerocoeff_select{3};            % Slope of CL-alpha curve [ ]

% Longitudinal stability
coef_fd = {-0.5347, -1.1494};
coef_ref = {-0.5718, -1.1935};
coef_select = coef_ref;

Cma    = coef_select{1};            % longitudinal stabilty [ ]
Cmde   = coef_select{2};            % elevator effectiveness [ ]

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

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.095;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
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
CZdt = 0.000;   %compared to the actual eleV0ator is an order of 
Cmdt = 0.000;   %magnitude smaller


a_11 = (CYbdot-2*mub)*b/V0;
a_22 = -b/(2*V0);
a_33 = -2*mub*KX2*(b/V0)^2;
a_34 = 2*mub*KXZ*(b/V0)^2;
a_41 = Cnbdot*b/V0;
a_43 = 2*mub*KXZ*(b/V0)^2;
a_44 = -2*mub*KZ2*(b/V0)^2;




% sys=ss(A_s,B_s(:,1),C,D(:,1));
% transf = tf(sys);
% e_A_s = eig(A_s);
% e_A_a = eig(A_a);
% time = linspace(0,145,1451);
% plot(time,p(index:index+1450));

