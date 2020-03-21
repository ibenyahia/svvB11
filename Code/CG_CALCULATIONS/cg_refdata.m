%NOTE: everything right now is in kg and m!!!
M_fuelused = 0.453592*([[360], 
    [412],
    [447],
    [478],
    [532],
    [570],
    [664],
    [694],
    [730],
    [755],
    [798],
    [825],
    [846],
    [881],
    [910]]);
    %Total fuel used

t = ([1157,1297,1426,1564,1787,1920,2239,2351,2484,2576,2741,2840,2920,3062,3166]);
%Measured time instances

M_BEM = 9165*0.453592*ones(15,1);               %Basic empty  mass
x_cgBEM = 291.65*2.54/100;                      %Cg position of the empty weight
                                                %given in appendix E?? NOT
                                                %SURE ABOUT THIS ONE
                                                
momentBEM = (M_BEM*x_cgBEM).*ones(15,1);

%Calculations for Payload cg
%Before seat change:
M_pilot1 = 95;
M_pilot2 = 92;
M_coordinator = 74;
M_1L = 66;
M_1R = 61;
M_2L = 75;
M_2R = 78;
M_3L = 86;
M_3R = 68;
M_Payload = (M_pilot1+M_pilot2+M_coordinator+M_1L+M_1R+M_2L+M_2R+M_3L+M_3R)*ones(15,1);
x_cgpilots = 131*2.54/100;
x_cg1 = 214*2.54/100;
x_cg2 = 251*2.54/100;
x_cg3 = 288*2.54/100;
x_cgcoordinator = 170*2.54/100;

momentPilots = ((M_pilot1+M_pilot2)*x_cgpilots);
momentCoordinator = M_coordinator*x_cgcoordinator;
moment1 = (M_1L+M_1R)*x_cg1;
moment2 = (M_2L+M_2R)*x_cg2;
moment3 = (M_3L+M_3R)*x_cg3;

momentPayload = (momentPilots+momentCoordinator+moment1+moment2+moment3)*ones(14,1);

%After seat 3R change:
x_cg3L = 288*2.54/100;   %not stated if either 3L of 3R moves
x_cg3R = 134*2.54/100;
moment3L = M_3L*x_cg3L;
moment3R = M_3R*x_cg3R;

momentPayload(15) = (momentPilots+momentCoordinator+moment1+moment2+moment3L+moment3R);


%Calculation for Fuel cg
m_i_fuel = 4050*0.453592;               %Initial fuel mass
M_i_fuel = ones(15,1)*m_i_fuel;               
M_fuel = M_i_fuel-M_fuelused;
xcg_fuel = 285.3*2.54/100;              %meters and assumed constant!
                                        %Because we (and the ref data) have
                                        %the first flight of the day.
                                        %Therefore, looking at table E2,
                                        %you can calculate the moment arm
                                        %at each mass given and see that
                                        %the changes between 3000-4100 lbs
                                        %are soooo small i think it can be
                                        %assumed constant. (285.135-285.5).
                                        %Averaged it equals = 285.3
momentFuel = M_fuel*xcg_fuel;


%Calulation of the total cg
x_cg = (momentBEM + momentPayload + momentFuel)./(M_BEM + M_Payload + M_fuel);

plot(x_cg, t-t(1));

