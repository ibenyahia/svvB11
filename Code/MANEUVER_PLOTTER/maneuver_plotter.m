clc
close all
%clear all
load('FLIGHT_and_REFdata\FTISxprt-20200309_flight1.mat')
% load('refdata.mat')
%get data from FTI file
t = flightdata.time.data ;                       %time [sec]
alpha = flightdata.vane_AOA.data ;               %angle of attack [deg]
dalpha_dt =flightdata.Ahrs1_bPitchRate.data;     % pitch rate [deg/s]
pitch = flightdata.Dadc1_bcAlt.data ;            % ptch [deg]
delta_e = flightdata.delta_e.data ;              % delt_e [deg] 


%Short Period
   %find begin time
h_sp =1 ;                                        %hours time stemp [h]
min_sp = 0 ;                                    %minutes time stamp [min]
sec_sp = 35 ;                                    %seconds time stamp [sec]
t_st_sp = 60^2 * h_sp + 60 * min_sp + sec_sp ;   % start time stamp [sec]

duar_min_sp = 0 ;                                % duaration maneuver [min]
duar_sec_sp = 10 ;                               % duaration maneuver [sec]
duar_tot_sp = duar_min_sp *60 + duar_sec_sp ;    %total duaration man [sec]
t_end_sp = t_st_sp + duar_tot_sp ;               % end time stamp [sec]

    %determine index
[dif_sp_0,i_sp_0] = min(abs(t-t_st_sp)) ;        %gives index of start point
[dif_sp_1,i_sp_1] = min(abs(t-t_end_sp)) ;       %gives index of end point

figure(1);
subplot(2,2,1);
plot(t(i_sp_0:i_sp_1)-t(i_sp_0),alpha(i_sp_0:i_sp_1))
title('alpha vs t')
ylabel('[deg]')
xlabel('[sec]')

subplot(2,2,3);
plot(t(i_sp_0:i_sp_1)-t(i_sp_0),dalpha_dt(i_sp_0:i_sp_1))
title('pitch rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')

subplot(2,2,2);
plot(t(i_sp_0:i_sp_1)-t(i_sp_0),delta_e(i_sp_0:i_sp_1))
title('delta_e rate vs t')
ylabel('[deg]')
xlabel('[sec]')

subplot(2,2,4);
plot(t(i_sp_0:i_sp_1)-t(i_sp_0),pitch(i_sp_0:i_sp_1))
title('pitch vs t')
ylabel('[deg]')
xlabel('[sec]')



%Phugoid
    %find begin time
h_ph =0 ;                                        %hours time stemp [h]
min_ph = 53 ;                                    %minutes time stamp [min] 
sec_ph = 57 ;                                    %seconds time stamp [sec]
t_st_ph= 60^2 * h_ph + 60 * min_ph + sec_ph ;    %start time stamp [sec]

duar_min_ph = 2 ;                                % duaration maneuver [min]
duar_sec_ph = 30 ;                               % duaration maneuver [sec]
duar_tot_ph = duar_min_ph *60 + duar_sec_ph ;    %total duaration man [sec]
t_end_ph = t_st_ph + duar_tot_ph ;               % end time stamp [sec]

    %determine index
[dif_ph_0,i_ph_0] = min(abs(t-t_st_ph)) ;        %gives index of start point
[dif_ph_1,i_ph_1] = min(abs(t-t_end_ph)) ;       %gives index of end point

figure(2);
subplot(2,1,1);
plot(t(i_ph_0:i_ph_1)-t(i_ph_0),alpha(i_ph_0:i_ph_1))
title('alpha vs t')
ylabel('[deg]')
xlabel('[sec]')

subplot(2,1,2);
plot(t(i_ph_0:i_ph_1)-t(i_ph_0),dalpha_dt(i_ph_0:i_ph_1))
title('ptch rate vs t')
ylabel('[deg/sec]')
xlabel('[sec]')

%Aperiodic Roll
    %find begin time
h_ar =0 ;                                       %hours time stemp [h]
min_ar = 59 ;                                   %minutes time stamp [min]
sec_ar = 10 ;                                   %seconds time stamp [sec]
t_st_ar = 60^2 * h_ar + 60 * min_ar + sec_ar ;  %start time stamp [sec]

duar_min_ar = 2 ;                               % duaration maneuver [min]
duar_sec_ar = 30 ;                              % duaration maneuver [sec]
duar_tot_ar = duar_min_ar *60 + duar_sec_ar ;   %total duaration man [sec]
t_end_ar = t_st_ar + duar_tot_ar ;              % end time stamp [sec]

    %determine index
[dif_ar_0,i_ar_0] = min(abs(t-t_st_ar)) ;        %gives index of start point
[dif_ar_1,i_ar_1] = min(abs(t-t_end_ar)) ;       %gives index of end point



%Spriral
    %find begin time
h_s =1 ;                                        %hours time stemp [h]
min_s = 5 ;                                    %minutes time stamp [min]
sec_s = 20 ;                                    %seconds time stamp [sec]
t_st_s = 60^2 * h_s + 60 * min_s + sec_s ;     %start time stamp [sec]

duar_min_s = 2 ;                                % duaration maneuver [min]
duar_sec_s = 30 ;                               % duaration maneuver [sec]
duar_tot_s = duar_min_s *60 + duar_sec_s ;      %total duaration man [sec]
t_end_s = t_st_s + duar_tot_s ;                 % end time stamp [sec]

    %determine index
[dif_s_0,i_s_0] = min(abs(t-t_st_s)) ;          %gives index of start point
[dif_s_1,i_s_1] = min(abs(t-t_end_s)) ;         %gives index of end point



%Dutch Roll
    %find begin time
h_dr =1 ;                                       %hours time stemp [h]
min_dr = 1 ;                                   %minutes time stamp [min]
sec_dr = 57 ;                                   %seconds time stamp [sec]
t_st_dr = 60^2 * h_dr + 60 * min_dr + sec_dr ; %start time stamp [sec]

duar_min_dr = 2 ;                               % duaration maneuver [min]
duar_sec_dr = 30 ;                              % duaration maneuver [sec]
duar_tot_dr = duar_min_dr *60 + duar_sec_dr ;   %total duaration man [sec]
t_end_dr = t_st_dr + duar_tot_dr ;              % end time stamp [sec]

    %determine index
[dif_dr_0,i_dr_0] = min(abs(t-t_st_dr)) ;        %gives index of start point
[dif_dr_1,i_dr_1] = min(abs(t-t_end_dr)) ;       %gives index of end point

%Dutch Roll damped
    %find begin time
h_dr_dp =1 ;                                       %hours time stemp [h]
min_dr_dp = 2 ;                                   %minutes time stamp [min]
sec_dr_dp = 47 ;                                   %seconds time stamp [sec]
t_st_dr_dp = 60^2 * h_dr_dp + 60 * min_dr_dp + sec_dr_dp ; %start time stamp [sec]

duar_min_dr_dp = 2 ;                               % duaration maneuver [min]
duar_sec_dr_dp = 30 ;                              % duaration maneuver [sec]
duar_tot_dr_dp = duar_min_dr_dp *60 + duar_sec_dr_dp ;   %total duaration man [sec]
t_end_dr_dp = t_st_dr_dp + duar_tot_dr_dp ;              % end time stamp [sec]
    
    %determine index
[dif_dr_dp_0,i_dr_dp_0] = min(abs(t-t_st_dr_dp)) ;       %gives index of start point
[dif_dr_dp_1,i_dr_dp_1] = min(abs(t-t_end_dr_dp)) ;      %gives index of end point




%a = flightdata.vane_AOA.data;
%t = flightdata.time.data;
%plot(t,a)