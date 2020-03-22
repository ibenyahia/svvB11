load('FTISxprt-20200309_flight1.mat')
i_SPfd = 30401;
t =  flightdata.time.data;
t_0 =t(30401);

[i_0,i_1]= indices(0,0,t_0,10,t);

theta = flightdata.Ahrs1_Pitch.data(i_0:i_1);
theta_dt = flightdata.Ahrs1_bPitchRate.data(i_0:i_1);
delta = flightdata.delta_e.data(i_0:i_1);


short_period_plot(theta,delta,theta_dt,t(i_0:i_1)-t_0,0,0,0,0);

