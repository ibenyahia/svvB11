%Short Period
   %find begin time
h_sp =0 ;                                        %hours time stemp [h]
min_sp = 50 ;                                    %minutes time stamp [min]
sec_sp = 50 ;                                    %seconds time stamp [sec]
t_tot_sp = 60^2 * h_sp + 60 * min_sp + sec_sp ;  % total time stamp [sec]

duar_min_sp = 2 ;                                % duaration maneuver [min]
duar_sec_sp = 30 ;                               % duaration maneuver [sec]
duar_tot_sp = duar_min_sp *60 + duar_sec_sp ;    %total duaration man [sec]
%Phugoid
    %find begin time
h_ph =0 ;                                        %hours time stemp [h]
min_ph = 54 ;                                    %minutes time stamp [min] 
sec_ph = 20 ;                                    %seconds time stamp [sec]
t_tot_ph= 60^2 * h_ph + 60 * min_ph + sec_ph ;   %total time stamp [sec]

duar_min_ph = 2 ;                                % duaration maneuver [min]
duar_sec_ph = 30 ;                               % duaration maneuver [sec]
duar_tot_ph = duar_min_ph *60 + duar_sec_ph ;     %total duaration man [sec]


%Aperiodic Roll
    %find begin time
h_ar =0 ;                                       %hours time stemp [h]
min_ar = 52 ;                                   %minutes time stamp [min]
sec_ar = 50 ;                                   %seconds time stamp [sec]
t_tot_ar = 60^2 * h_ar + 60 * min_ar + sec_ar ; %total time stamp [sec]

duar_min_ar = 2 ;                               % duaration maneuver [min]
duar_sec_ar = 30 ;                              % duaration maneuver [sec]
duar_tot_ar = duar_min_ar *60 + duar_sec_ar ;    %total duaration man [sec]


%Spriral
    %find begin time
h_s =1 ;                                        %hours time stemp [h]
min_s = 1 ;                                    %minutes time stamp [min]
sec_s = 0 ;                                    %seconds time stamp [sec]
t_tot_s = 60^2 * h_s + 60 * min_s + sec_s ;     %total time stamp [sec]

duar_min_s = 2 ;                                % duaration maneuver [min]
duar_sec_s = 30 ;                               % duaration maneuver [sec]
duar_tot_s = duar_min_s *60 + duar_sec_s ;      %total duaration man [sec]


%Dutch Roll
    %find begin time
h_dr =0 ;                                       %hours time stemp [h]
min_dr = 57 ;                                   %minutes time stamp [min]
sec_dr = 30 ;                                   %seconds time stamp [sec]
t_tot_dr = 60^2 * h_dr + 60 * min_dr + sec_dr ; %total time stamp [sec]

duar_min_dr = 2 ;                               % duaration maneuver [min]
duar_sec_dr = 30 ;                              % duaration maneuver [sec]
duar_tot_dr = duar_min_dr *60 + duar_sec_dr ;   %total duaration man [sec]

%Dutch Roll damped
    %find begin time
h_dr_dp =0 ;                                       %hours time stemp [h]
min_dr_dp = 58 ;                                   %minutes time stamp [min]
sec_dr_dp = 50 ;                                   %seconds time stamp [sec]
t_tot_dr_dp = 60^2 * h_dr_dp + 60 * min_dr_dp + sec_dr_dp ; %total time stamp [sec]

duar_min_dr_dp = 2 ;                               % duaration maneuver [min]
duar_sec_dr_dp = 30 ;                              % duaration maneuver [sec]
duar_tot_dr_dp = duar_min_dr_dp *60 + duar_sec_dr_dp ;   %total duaration man [sec]




a = flightdata.vane_AOA.data;
t = flightdata.time.data;
plot(t,a)