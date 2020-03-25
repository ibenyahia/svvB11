% Find mass at stationary measurement timesteps for flight data

M_fuelused = 0.453592*([[443], 
    [468],
    [527],
    [555],
    [577],
    [600],
    [686],
    [698],
    [711],
    [730],
    [757],
    [787],
    [815],
    [838],
    [865]]);
    %Total fuel used

t = ([1350,1424,1651,1765,1865,1966,2215,2280,2328,2415,2530,2672,2764,2860,2962]);
%Measured time instances

M_BEM = 9165*0.453592;                  %Basic empty  mass in kg
M_pax = 740;                            %Pax mass in kg
M_i_fuel = 4100*0.453592;               %Initial fuel mass
m_i_ramp = M_BEM+M_pax+M_i_fuel;        %Initial ramp mass
M_i_ramp = ones(15,1)*m_i_ramp;         %Matrix of initial ramp mass

Mtotal = M_i_ramp - M_fuelused;         %Total mass of the aircraft during flight

plot(t-t(1), Mtotal,'b')
xlabel('Time [s]')
ylabel('Total mass [kg]')
title('Total mass of the aircraft of reference flight')
