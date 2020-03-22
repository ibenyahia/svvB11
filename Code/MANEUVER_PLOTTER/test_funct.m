load('FTISxprt-20200309_flight1.mat')
i_DRfd = 34411;
t =  flightdata.time.data;
t_0 =t(i_DRfd);

[i_0,i_1]= indices(0,0,t_0,25,t);

p = flightdata.Ahrs1_bRollRate.data(i_0:i_1);      % roll rate [rad/s]
r = flightdata.Ahrs1_bYawRate.data(i_0:i_1);       % yaw rate [rad/s]
phi = flightdata.Ahrs1_Roll.data(i_0:i_1);         % roll [rad]
 

theta = flightdata.Ahrs1_Pitch.data(i_0:i_1);
theta_dt = flightdata.Ahrs1_bPitchRate.data(i_0:i_1);
delta = flightdata.delta_e.data(i_0:i_1);


dr_plot(p,phi,r,t(i_0:i_1)-t_0,0,0,0,0,0);

