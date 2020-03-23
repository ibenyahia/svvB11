load('FLIGHT_and_REFdata\FTISxprt-20200309_flight1.mat')

% Get variables 
t = flightdata.time.data;                        % time [sec]
hp = flightdata.Dadc1_alt.data*0.3048;           % altitude [m]
alpha = flightdata.vane_AOA.data*pi/180 ;        % angle of attack [rad]
theta = flightdata.Ahrs1_Pitch.data*pi/180 ;     % pitch [rad]
vtas = flightdata.Dadc1_tas.data*0.514444;       % true airspeed [m/s]  

% Index for flight data

i_SPfd = 30401;
i_PHfd = 32511;
i_DRfd = 34411;
i_DRDampfd = 35211;
i_APfd = 31611;
i_SPIRfd = 36781; % 36511

index = i_SPIRfd; % Index of interest for initial conditions

%Order in list: hp0 [m], V0[m/s], alpha0 [rad], th0 [rad], mass [kg]

%Flight data
shortperiod_fd = {hp(i_SPfd),vtas(i_SPfd),alpha(i_SPfd),theta(i_SPfd),Mtotal(i_SPfd)};
phugoid_fd = {hp(i_PHfd),vtas(i_PHfd),alpha(i_PHfd),theta(i_PHfd),Mtotal(i_PHfd)};
dr_fd = {hp(i_DRfd),vtas(i_DRfd),alpha(i_DRfd),theta(i_DRfd),Mtotal(i_DRfd)};
drdamp_fd = {hp(i_DRDampfd),vtas(i_DRDampfd),alpha(i_DRDampfd),theta(i_DRDampfd),Mtotal(i_DRDampfd)};
apr_fd = {hp(i_APfd),vtas(i_APfd),alpha(i_APfd),theta(i_APfd),Mtotal(i_APfd)};
spiral_fd = {hp(i_SPIRfd),vtas(i_SPIRfd),alpha(i_SPIRfd),theta(i_SPIRfd),Mtotal(i_SPIRfd)};

selection = spiral_fd;   %Replace name with flight condition of interest

% Initial conditions
hp0    = selection{1}(1);             % Initial height [m]
V0     = selection{2}(1);             % Initial airspeed [m/s]
alpha0 = selection{3}(1);                           % Alpha stability-axis [rad]
th0    = selection{4}(1);                           % Theta stability-axis [rad]
m      = selection{5}(1);         	  % Initial aircraft mass [kg]

e      = 0.725 ;            % Oswald factor [-]
CD0    = 0.0214 ;            % Zero lift drag coefficient [-]
CLa    = 4.59 ;            % Slope of CL-alpha curve [rad^-1]

% Longitudinal stability
Cma    = -0.5347;            % longitudinal stabilty [rad^-1]
Cmde   = -1.1494;            % elevator effectiveness [rad^-1]

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

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)

W      = m*g; % [N]       (aircraft weight)

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

CL    = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]

CD    = CD0 + CL^2/(pi*A*e);  % Drag coefficient [ ]

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

%Symmetric eigenmotions

%Short period motion
A_SPM = 2*muc*KY2*(2*muc-CZadot);
B_SPM = -2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq;
C_SPM = CZa*Cmq-(2*muc+CZq)*Cma;

Lambda_1SPM = (-B_SPM+sqrt(-4*A_SPM*C_SPM+B_SPM^2))/(2*A_SPM)*V0/c;
Lambda_2SPM = (-B_SPM-sqrt(-4*A_SPM*C_SPM+B_SPM^2))/(2*A_SPM)*V0/c;
Lambda_scaSPM = sqrt(real(Lambda_1SPM)^2+imag(Lambda_1SPM)^2);

P_SPM = (2*pi/imag(Lambda_1SPM))*(c/V0);
T_SPM = (log(1/2)/Lambda_scaSPM)*(c/V0);
D_SPM = -real(Lambda_1SPM)/(Lambda_scaSPM);

%Phugoid
A_PHU = 2*muc*(CZa*Cmq-2*muc*Cma);
B_PHU = 2*muc*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa);
C_PHU = CZ0*(Cmu*CZa-CZu*Cma);

Lambda_1PH = (-B_PHU+sqrt(-4*A_PHU*C_PHU+B_PHU^2))/(2*A_PHU)*V0/c;
Lambda_2PH = (-B_PHU-sqrt(-4*A_PHU*C_PHU+B_PHU^2))/(2*A_PHU)*V0/c;
Lambda_scaPH = sqrt(real(Lambda_1PH)^2+imag(Lambda_1PH)^2);

P_PHU = (2*pi/imag(Lambda_1PH))*(c/V0);
T_PHU = (log(1/2)/Lambda_scaPH)*(c/V0);
D_PHU = -real(Lambda_1PH)/(Lambda_scaPH);


%%

%Asymmetric eigenmotions

%Aperiodic roll motion

Lambda_AP = Clp/(4*mub*KX2)*V0/b;
Lambda_scaAP = sqrt(real(Lambda_AP)^2+imag(Lambda_AP)^2);

P_AP = (2*pi/imag(Lambda_AP))*(b/V0);
T_AP = (log(1/2)/Lambda_scaAP)*(b/V0);
D_AP = -real(Lambda_AP)/(Lambda_scaAP);


%Dutch roll

Lambda_1DR = 2*(Cnr+2*KZ2*CYb)+sqrt(-64*KZ2*(4*mub*Cnb+CYb*Cnr)+4*(Cnr+2*KZ2*CYb)^2)/(16*mub*KZ2)*V0/b;
Lambda_2DR = 2*(Cnr+2*KZ2*CYb)-sqrt(-64*KZ2*(4*mub*Cnb+CYb*Cnr)+4*(Cnr+2*KZ2*CYb)^2)/(16*mub*KZ2)*V0/b;
Lambda_scaDR = sqrt(real(Lambda_1DR)^2+imag(Lambda_1DR)^2);

P_DR = (2*pi/imag(Lambda_1DR))*(b/V0);
T_DR = (log(1/2)/Lambda_scaDR)*(b/V0);
D_DR = -real(Lambda_1DR)/(Lambda_scaDR);

%Spiral motion

Lambda_SPIR = (2*CL*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))*V0/b;
Lambda_scaSPIR = sqrt(real(Lambda_SPIR)^2+imag(Lambda_SPIR)^2);

P_SPIR = (2*pi/imag(Lambda_SPIR))*(b/V0);
T_SPIR = (log(1/2)/Lambda_scaSPIR)*(b/V0);
D_SPIR = -real(Lambda_SPIR)/(Lambda_scaSPIR);

