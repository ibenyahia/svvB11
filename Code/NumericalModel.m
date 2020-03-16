% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

hp0    = 3000;      	  % pressure altitude in the stationary flight condition [m]
V0     = 72;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = 0.0872665;       	  % angle of attack in the stationary flight condition [rad]
th0    = 0.0349066;       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 3500;         	  % mass [kg]

% aerodynamic properties
e      = 0.8;            % Oswald factor [ ]
CD0    = 0.04;            % Zero lift drag coefficient [ ]
CLa    = 5.084;            % Slope of CL-alpha curV0e [ ]

% Longitudinal stability
Cma    = -0.5626;            % longitudinal stabilty [ ]
Cmde   = -1.1642;            % elevator effectiveness [ ]

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

% Constant V0alues concerning aircraft inertia

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
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity deriV0atiV0es

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





%Numerical model stuff starts here
%V0ARIABLES THAT ARE MISSING


CXdt = 0.000;   %trim tab is really small for this aircraft so it's contribution
CZdt = 0.000;   %compared to the actual eleV0ator is an order of 
Cmdt = 0.000;   %magnitude smaller


x_u = (V0/c) * CXu/(2*muc);
x_alpha = (V0/c) * CXa/(2*muc);
x_theta = (V0/c) * CZ0/(2*muc);
x_q = (V0/c) * CXq/(2*muc);
x_deltae = (V0/c) * CXde/(2*muc);
x_deltat = (V0/c) * CXdt/(2*muc);

z_u = (V0/c) * CZu/(2*muc - CZadot);
z_alpha = (V0/c) * CZa/(2*muc - CZadot);
z_theta = - (V0/c) * CX0/(2*muc - CZadot);
z_q = (V0/c) * (2*muc+CZq)/(2*muc - CZadot);
z_deltae = (V0/c) * CZde/(2*muc - CZadot);   
z_deltat = (V0/c) * CZdt/(2*muc - CZadot);

m_u = (V0/c) * (Cmu + CZu * (Cmadot/(2*muc-CZadot)))/(2*muc*KY2);
m_alpha = (V0/c) * (Cma + CZa * (Cmadot/(2*muc-CZadot)))/(2*muc*KY2);
m_theta = -(V0/c) * (CX0 * (Cmadot/(2*muc-CZadot)))/(2*muc*KY2);
m_q = (V0/c) * (Cmq + Cmadot * ((2*muc+CZq)/(2*muc-CZadot)))/(2*muc*KY2);
m_deltae = (V0/c) * (Cmde + CZde * (Cmadot/(2*muc-CZadot)))/(2*muc*KY2); 
m_deltat = (V0/c) * (Cmdt + CZdt * (Cmadot/(2*muc-CZadot)))/(2*muc*KY2);

y_beta = (V0/b) * CYb/(2*mub);
y_phi = (V0/b) * CL/(2*mub); 
y_p = (V0/b) * CYp/(2*mub);
y_r = (V0/b) * (CYr-4*mub)/(2*mub);
y_deltaa = (V0/b) * CYda/(2*mub);
y_deltar = (V0/b) * CYdr/(2*mub);

l_beta = (V0/b) * (Clb*KZ2+Cnb*KXZ)/(4*mub*(KX2*KZ2-KXZ^2));
l_phi = 0;
l_p = (V0/b) * (Clp*KZ2+Cnp*KXZ)/(4*mub*(KX2*KZ2-KXZ^2)); 
l_r = (V0/b) * (Clr*KZ2+Cnr*KXZ)/(4*mub*(KX2*KZ2-KXZ^2));
l_deltaa = (V0/b) * (Clda*KZ2+Cnda*KXZ)/(4*mub*(KX2*KZ2-KXZ^2));
l_deltar = (V0/b) * (Cldr*KZ2+Cndr*KXZ)/(4*mub*(KX2*KZ2-KXZ^2));

n_beta = (V0/b) * (Clb*KXZ+Cnb*KX2)/(4*mub*(KX2*KZ2-KXZ^2));
n_phi = 0;
n_p = (V0/b) * (Clp*KXZ+Cnp*KX2)/(4*mub*(KX2*KZ2-KXZ^2));
n_r = (V0/b) * (Clr*KXZ+Cnr*KX2)/(4*mub*(KX2*KZ2-KXZ^2));
n_deltaa = (V0/b) * (Clda*KXZ+Cnda*KX2)/(4*mub*(KX2*KZ2-KXZ^2));
n_deltar = (V0/b) * (Cldr*KXZ+Cndr*KX2)/(4*mub*(KX2*KZ2-KXZ^2));


A_s = [x_u, x_alpha, x_theta, 0;...
    z_u, z_alpha, z_theta, z_q;...
    0, 0, 0, V0/c;...
    m_u, m_alpha, m_theta, m_q];

B_s = [x_deltae, x_deltat;...
    z_deltae, z_deltat;...
    0, 0;...
    m_deltae, m_deltat];


A_a = [y_beta, y_phi, y_p, y_r;...
    0, 0, 2*V0/b, 0;...
    l_beta, 0, l_p, l_r;...
    n_beta, 0, n_p, n_r;];

B_a = [0, y_deltar;...
    0, 0;...
   l_deltaa, l_deltar;...
   n_deltaa, n_deltar];

C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];

D = [0,0;...
    0,0;...
    0,0;...
    0,0];

sys=ss(A_a,B_a,C,D);

e_A_s = eig(A_s);
e_A_a = eig(A_a);

subplot(2,1,1)
step(sys)
subplot(2,1,2)
impulse(sys)