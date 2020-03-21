% V = 0.514444444*([[249],    %in m/s but IAS and not EAS
%     [221],
%     [192],
%     [163],
%     [130],
%     [118],
%     [161],
%     [150],
%     [140],
%     [130],
%     [173],
%     [179],
%     [192],
%     [161],
%     [161]]);
% 
% H = 0.3048*([[5010],                %Height in meters
%     [5020],
%     [5020],
%     [5030],
%     [5020],
%     [5110],
%     [6060],
%     [6350],
%     [6550],
%     [6880],
%     [6160],
%     [5810],
%     [5310],
%     [5730],
%     [5790]]);
% 
% TAT = 273.15 + ([[12.5],
%     [10.5],
%     [8.8],
%     [7.2],
%     [6],
%     [5.2],
%     [5.5],
%     [4.5],
%     [3.5],
%     [2.5],
%     [5],
%     [6.2],
%     [8.2],
%     [5],
%     [5]]);
% 
% rho = ([[1.05523],                 %individual ISA calculations performed
%     [1.05491],
%     [1.05491],
%     [1.05459],
%     [1.05491],
%     [1.05203],
%     [1.02205],
%     [1.01304],
%     [1.00685],
%     [0.996704],
%     [1.01894],
%     [1.02988],
%     [1.04567],
%     [1.03239],
%     [1.03051]]);                  
% 
% M_fuelused = ([[163.29], %Fuel used in kg
%     [186.88],
%     [202.76],
%     [216.82],
%     [241.31],
%     [258.55],
%     [301.19],
%     [314.70],
%     [331.13],
%     [342.46],
%     [361.97],
%     [374.21],
%     [383.74],
%     [399.61],
%     [412.77]]);
% 
% t = ([1157,1297,1426,1564,1787,1920,2239,2351,2484,2576,2741,2840,2920,3062,3166]);
% %Measured time instances
% 
% S =  30.00 * ones(15,1);                %The wing area
% 
% M_BEM = 9165*0.453592;                  %Basic empty  mass in kg
% M_pax = 695;                            %Pax mass in kg
% M_i_fuel = 4050*0.453592;               %Initial fuel mass
% m_i_ramp = M_BEM+M_pax+M_i_fuel;        %Initial ramp mass
% M_i_ramp = ones(15,1)*m_i_ramp;         %Matrix of initial ramp mass
% 
% Mtotal = M_i_ramp - M_fuelused;         %Total mass of the aircraft during flight
% 
% CL_config1 = 9.81.*Mtotal(1:6)./(1/2.*rho(1:6).*V(1:6).^2.*S(1:6));
% CL_config2 = 9.81.*Mtotal(7:15)./(1/2.*rho(7:15).*V(7:15).^2.*S(7:15));
% 
% a = ([[1.7],
%     [2.4],
%     [3.6],
%     [5.4],
%     [8.7],
%     [10.6],
%     [5.3],
%     [6.3],
%     [7.3],
%     [8.5],
%     [4.5],
%     [4.1],
%     [3.4],
%     [5.3],
%     [5.3]]);
% 
% 
% %a1 = plot(a(1:6),CL_config1, 'o-r');
% %M1 = 'Configuration 1';
% %xlabel('Angle of attack [deg]')
% %ylabel('CL [-]')
% %title('CL-alpha curve from configuration 1 and 2')
% 
% %hold on
% 
% %a2 = plot(a(7:15), CL_config2, 'o-b');
% %M2 = 'Configuration 2';
% 
% %legend([a1;a2], M1, M2)
% %hold off
% 
% 
% %CLa for the first configuration
% %Ff de goede punten nog tweaken want vooral AOA gaat beetje door elkaar
% Cla_config1 = (CL_config1(6)-CL_config1(1))/(a(6)-a(1));
% 
% %CLa for the second configuration
% %Ff de goede punten nog tweaken want vooral AOA gaat beetje door elkaar
% CLa_config2 = (CL_config2(4)-CL_config2(2))/(a(10)-a(8));
% 
% 
% 
% de = ([[0],
%     [-0.4],
%     [-0.9],
%     [-1.5],
%     [0.4],
%     [0.6],
%     [1],
%     [0],
%     [-0.5]]);
% 
% %plot(a(7:15), de, 'o-b')
% %xlabel('Angle of attack [deg]')
% %ylabel('Elevator deflection [deg]')
% 
% slope_de_a = (de(7)-de(4))/(a(13)-a(10)); %ff goed kijken welke waardes om te 
%                                          %kiezen voor de slope calculation!
%                                          %ze gaan namelijk door elkaar heen
%                                          %bijv bij de elevator heel erg.
% 
% 
%                                          
%                                          
% %CMdelta calculation:
% c = 2.0569;        %Chord length
% cg_change = 7.1462-7.1036;        
% Cmd = -(1/(de(9)-de(8)))*(CL_config2(9))*(cg_change/c)
% 
% %CMalpha calculation
% CMa = -Cmd*slope_de_a




