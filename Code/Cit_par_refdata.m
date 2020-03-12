% Citation 550 - Linear simulation hai

% xcg = 0.25*c

% Stationary flig
% ht condition


hp0    = 0.3048*[5010,5020,5020,5030,5020,5110];      	  % pressure altitude in the stationary flight condition [m]
V0     = 0.514444*([249,221,192,163,130,118]-2);            % calibrated airspeed in the stationary flight condition [m/sec]
alpha0 = [1.7,2.4,3.6,5.4,8.7,10.6];       	  % angle of attack in the stationary flight condition [rad]
th0    = alpha0;       	  % pitch angle in the stationary flight condition [rad]
T_m = 273.15 + [12.5,10.5,8.8,7.2,6,5.2];
FFl = [798,673,561,463,443,474]/3600*0.45359237;
FFr = [813,682,579,484,467,499]/3600*0.45359237;

[rho,p,V_EAS,V_TAS,delta_T,M_T] = getImpValues(hp0,V0,T_m);
% Aircraft mass
m      = Mtotal;         	  % mass [kg]
%%

% aerodynamic properties
% e      = ;            % Oswald factor [ ]
% CD0    = ;            % Zero lift drag coefficient [ ]
% CLa    = ;            % Slope of CL-alpha curve [ ]
% 
% % Longitudinal stability
% Cma    = ;            % longitudinal stabilty [ ]
% Cmde   = ;            % elevator effectiveness [ ]

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
W = reshape(m*g,1,length(m));             % [N]       (aircraft weight)
Thr_L = [3176.3,2564.84,2006.48,1508.82,1549.65,1832.57];
Thr_R = [3271.58,2623.89,2127.61,1654.57,1720.84,2032.14];
Thr_tot = Thr_L + Thr_R;
CD = Thr_tot./(0.5*rho0*V_EAS.^2*S);
CL = W(1:6)./(0.5*rho0*V_EAS.^2*S);
CL_2 = CL.^2;
pAe = 23.8; % From curve fit
e = pAe/(pi*A);
CD0 = 0.425/(pAe);
plot(CD,CL_2);

%%
% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

% CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
% CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.02792;
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

As = array([])
