% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition
                    %(SPM, PHU, ARM, DR, DRDAMP, SP)
hp0    = 0.3048*([10227, 10070, 10042, 10107, 10087, 11808]);      	  % pressure altitude in the stationary flight condition [m]
V0     = 0.514444*([179.5590, 176.2233, 183.6735, 181.6761, 181.5635, 178.3813]);            % true airspeed in the stationary flight condition [m/sec]
alpha0 = (pi/180)*([6.4449, 5.9677, 7.4384, 5.8394, 5.6176, 5.9624]);       	  % angle of attack in the stationary flight condition [rad]
th0    = (pi/180)*([5.9732, 4.8115, 0.5733, 4.6057, 3.9284, 4.8012]);       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = ([6345.7, 6331.5, 6341.6, 6310.5, 6301.6, 6280.2]);  	  % mass [kg]

% aerodynamic properties
e      = 0.742;            % Oswald factor [-]
CD0    = 0.0209;            % Zero lift drag coefficient [-]
CL0    = 0.0833;           % CL at 0 aoa [-]
CLa    = 0.084;            % Slope of CL-alpha curve [rad^-1]

% Longitudinal stability
Cma    = -0.5347;            % longitudinal stabilty [deg^-1]
Cmde   = -1.1494;            % elevator effectiveness [rad&]

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

rho_SPM     = rho0*((1+(lambda*hp0(1)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho_PHU     = rho0*((1+(lambda*hp0(2)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho_ARM     = rho0*((1+(lambda*hp0(3)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho_DR      = rho0*((1+(lambda*hp0(4)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho_DRDAMP  = rho0*((1+(lambda*hp0(5)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho_SM      = rho0*((1+(lambda*hp0(6)/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)

W_SPM      = m(1)*g;				                        % [N]       (aircraft weight)
W_PHU      = m(2)*g;				                        % [N]       (aircraft weight)
W_ARM      = m(3)*g;				                        % [N]       (aircraft weight)
W_DR       = m(4)*g;				                        % [N]       (aircraft weight)
W_DRDAMP   = m(5)*g;				                        % [N]       (aircraft weight)
W_SM       = m(6)*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc_SPM     = m(1)/(rho_SPM*S*c);
muc_PHU     = m(2)/(rho_PHU*S*c);

mub_ARM     = m(3)/(rho_ARM*S*b);
mub_DR      = m(4)/(rho_DR*S*b);
mub_DRDAMP  = m(5)/(rho_DRDAMP*S*b);
mub_SM      = m(6)/(rho_SM*S*b);

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

CL_SPM      = 2*W_SPM/(rho_SPM*(V0(1)^2)*S);               % Lift coefficient [ ]
CL_PHU      = 2*W_PHU/(rho_PHU*(V0(2)^2)*S);               % Lift coefficient [ ]
CL_ARM      = 2*W_ARM/(rho_ARM*(V0(3)^2)*S);               % Lift coefficient [ ]
CL_DR       = 2*W_DR/(rho_DR*(V0(4)^2)*S);               % Lift coefficient [ ]
CL_DRDAMP   = 2*W_DRDAMP/(rho_DRDAMP*(V0(5)^2)*S);               % Lift coefficient [ ]
CL_SM       = 2*W_SM/(rho_SM*(V0(6)^2)*S);               % Lift coefficient [ ]

CD_SPM      = CD0+((CL_SPM)^2)/(pi*A*e);  % Drag coefficient [ ]
CD_PHU      = CD0+((CL_PHU)^2)/(pi*A*e);  % Drag coefficient [ ]
CD_ARM      = CD0+((CL_ARM)^2)/(pi*A*e);  % Drag coefficient [ ]
CD_DR       = CD0+((CL_DR)^2)/(pi*A*e);  % Drag coefficient [ ]
CD_DRDAMP   = CD0+((CL_DRDAMP)^2)/(pi*A*e);  % Drag coefficient [ ]
CD_SM       = CD0+((CL_SM)^2)/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0_SPM     = W_SPM*sin(th0(1))/(0.5*rho_SPM*(V0(1)^2)*S);
CX0_PHU     = W_PHU*sin(th0(2))/(0.5*rho_PHU*(V0(2)^2)*S);
CX0_ARM     = W_ARM*sin(th0(3))/(0.5*rho_ARM*(V0(3)^2)*S);
CX0_DR      = W_DR*sin(th0(4))/(0.5*rho_DR*(V0(4)^2)*S);
CX0_DRDAMP  = W_DRDAMP*sin(th0(5))/(0.5*rho_DRDAMP*(V0(5)^2)*S);
CX0_SM      = W_SM*sin(th0(6))/(0.5*rho_SM*(V0(6)^2)*S);
CXu    = -0.095;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0_SPM     = -W_SPM*cos(th0(1))/(0.5*rho_SPM*(V0(1)^2)*S);
CZ0_PHU     = -W_PHU*cos(th0(2))/(0.5*rho_PHU*(V0(2)^2)*S);
CZ0_ARM     = -W_ARM*cos(th0(3))/(0.5*rho_ARM*(V0(3)^2)*S);
CZ0_DR      = -W_DR*cos(th0(4))/(0.5*rho_DR*(V0(4)^2)*S);
CZ0_DRDAMP  = -W_DRDAMP*cos(th0(5))/(0.5*rho_DRDAMP*(V0(5)^2)*S);
CZ0_SM      = -W_SM*cos(th0(6))/(0.5*rho_SM*(V0(6)^2)*S);
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
