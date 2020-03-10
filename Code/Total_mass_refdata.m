M_fuelused = ([[163.29], %Fuel used
    [186.88],
    [202.76],
    [216.82],
    [241.31],
    [258.55],
    [301.19],
    [314.70],
    [331.13],
    [342.46],
    [361.97],
    [374.21],
    [383.74],
    [399.61],
    [412.77]])

t = ([1157,1297,1426,1564,1787,1920,2239,2351,2484,2576,2741,2840,2920,3062,3166])
%Measured time instances

M_BEM = 9165*0.453592;                  %Basic empty  mass in kg
M_pax = 695;                            %Pax mass in kg
M_i_fuel = 4050*0.453592;               %Initial fuel mass
m_i_ramp = M_BEM+M_pax+M_i_fuel;        %Initial ramp mass
M_i_ramp = ones(15,1)*m_i_ramp;         %Matrix of initial ramp mass

Mtotal = M_i_ramp - M_fuelused;         %Total mass of the aircraft during flight

plot(t, Mtotal,'b')
xlabel('Time [s]')
ylabel('Total mass [kg]')
title('Total mass of the aircraft of reference flight')
