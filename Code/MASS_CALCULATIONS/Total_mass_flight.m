load('FLIGHT_and_REFdata\FTISxprt-20200309_flight1.mat')

m_le = flightdata.lh_engine_FU.data;    %Extract mass left engine in lbs    
m_re = flightdata.rh_engine_FU.data;    %Extract mass right engine in lbs
t = flightdata.time.data;               %Extract time of flight

M_BEM = 9165*0.453592;                  %Basic empty  mass in kg
M_pax = 740;                            %Pax mass in kg
M_i_fuel = 4100*0.453592;               %Initial fuel mass EXAMPLE MASS!
m_i_ramp = M_BEM+M_pax+M_i_fuel;        %Initial ramp mass
M_i_ramp = ones(49741,1)*m_i_ramp;      %Matrix of initial ramp mass

m_LE = m_le*0.453592;                   %Fuel used left engine during flight in kg
m_RE = m_re*0.453592;                   %Fuel used right engine during flight in kg
M_fuelused = m_LE+m_RE;                     %Total fuel used during flight in kg

Mtotal = M_i_ramp - M_fuelused;             %Total mass of the aircraft during flight

plot(t, Mtotal,'b')
xlabel('Time [s]')  
ylabel('Total mass [kg]')
% title('Total Mass of the Aircraft during Flight')
